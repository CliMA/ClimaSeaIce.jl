#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/../.." && pwd)"

CICE_ROOT="${CICE_ROOT:-/private/tmp/cice-validation/CICE}"
CASE_ROOT="${CASE_ROOT:-/private/tmp/cice-validation/cases}"
MACHINE="${CICE_MACHINE:-conda}"
ENVIRONMENT="${CICE_ENVIRONMENT:-homebrew}"
PES="${CICE_PES:-1x1}"
EXECUTE="${RUN_CICE_CASES_EXECUTE:-0}"
BUILD="${RUN_CICE_CASES_BUILD:-0}"
RUN="${RUN_CICE_CASES_RUN:-0}"
INSTALL_TEMPLATES="${CICE_INSTALL_LOCAL_TEMPLATES:-1}"
ALLOW_TEMPLATE_OVERWRITE="${CICE_ALLOW_TEMPLATE_OVERWRITE:-0}"
CASE_FILTER="${CICE_CASE_FILTER:-}"
CONDUCT_FILTER="${CICE_CONDUCT_FILTER:-}"
CASE_SUFFIX="${CICE_CASE_SUFFIX:-}"
INTERNAL_SNOW_THICKNESS="${CICE_INTERNAL_SNOW_THICKNESS:-}"
PRESCRIBED_WARMING_SURFACE_FLUX="${CICE_WARMING_SURFACE_FLUX_W_M2:-8.0}"
PRESCRIBED_ABLATION_SURFACE_FLUX="${CICE_ABLATION_SURFACE_FLUX_W_M2:-75.0}"
BALANCE_COLD_BOTTOM_FLUX="${CICE_BALANCE_COLD_BOTTOM_FLUX:-0}"

cases=(
  cold_conductive_relaxation
  surface_warming
  surface_ablation
  basal_growth
)

conducts=(MU71 bubbly)

runtime_set() {
  local case_id="$1"

  case "${case_id}" in
    surface_warming)
      echo run2day
      ;;
    *)
      echo run5day
      ;;
  esac
}

conduct_set() {
  local conduct="$1"

  case "${conduct}" in
    MU71)
      echo bl99_conduct_mu71
      ;;
    bubbly)
      echo bl99_conduct_bubbly
      ;;
    *)
      echo "Unknown conduct variant: ${conduct}" >&2
      return 2
      ;;
  esac
}

case_set() {
  local case_id="$1"

  case "${case_id}" in
    cold_conductive_relaxation|surface_warming|surface_ablation|basal_growth)
      echo "bl99_case_${case_id}"
      ;;
    *)
      echo "Unknown validation case: ${case_id}" >&2
      return 2
      ;;
  esac
}

install_template() {
  local source_file="$1"
  local target_file="$2"
  local target_dir

  target_dir="$(dirname "${target_file}")"
  mkdir -p "${target_dir}"

  if [[ ! -e "${target_file}" ]]; then
    cp "${source_file}" "${target_file}"
    echo "Installed ${target_file}"
  elif cmp -s "${source_file}" "${target_file}"; then
    echo "Already installed ${target_file}"
  elif [[ "${ALLOW_TEMPLATE_OVERWRITE}" == "1" ]]; then
    cp "${source_file}" "${target_file}"
    echo "Updated ${target_file}"
  else
    echo "Refusing to overwrite differing CICE file: ${target_file}" >&2
    echo "Set CICE_ALLOW_TEMPLATE_OVERWRITE=1 if this local checkout should use the validation template." >&2
    return 2
  fi
}

install_local_templates() {
  local template_root="${SCRIPT_DIR}/cice_templates"

  install_template "${template_root}/machines/env.conda_homebrew" \
                   "${CICE_ROOT}/configuration/scripts/machines/env.conda_homebrew"
  install_template "${template_root}/machines/Macros.conda_homebrew" \
                   "${CICE_ROOT}/configuration/scripts/machines/Macros.conda_homebrew"
  install_template "${template_root}/options/set_nml.col" \
                   "${CICE_ROOT}/configuration/scripts/options/set_nml.col"
  install_template "${template_root}/options/set_nml.bl99_validation" \
                   "${CICE_ROOT}/configuration/scripts/options/set_nml.bl99_validation"
  install_template "${template_root}/options/set_nml.bl99_conduct_mu71" \
                   "${CICE_ROOT}/configuration/scripts/options/set_nml.bl99_conduct_mu71"
  install_template "${template_root}/options/set_nml.bl99_conduct_bubbly" \
                   "${CICE_ROOT}/configuration/scripts/options/set_nml.bl99_conduct_bubbly"
  install_template "${template_root}/options/set_nml.bl99_case_cold_conductive_relaxation" \
                   "${CICE_ROOT}/configuration/scripts/options/set_nml.bl99_case_cold_conductive_relaxation"
  install_template "${template_root}/options/set_nml.bl99_case_surface_warming" \
                   "${CICE_ROOT}/configuration/scripts/options/set_nml.bl99_case_surface_warming"
  install_template "${template_root}/options/set_nml.bl99_case_surface_ablation" \
                   "${CICE_ROOT}/configuration/scripts/options/set_nml.bl99_case_surface_ablation"
  install_template "${template_root}/options/set_nml.bl99_case_basal_growth" \
                   "${CICE_ROOT}/configuration/scripts/options/set_nml.bl99_case_basal_growth"
}

patch_internal_snow_thickness() {
  local thickness="$1"
  local source_file="${CICE_ROOT}/cicecore/cicedyn/general/ice_init.F90"

  if [[ ! -f "${source_file}" ]]; then
    echo "CICE source file not found: ${source_file}" >&2
    return 2
  fi

  if ! grep -Eq "hsno_init = [0-9.]+_dbl_kind" "${source_file}"; then
    echo "Could not find hsno_init assignment in ${source_file}" >&2
    return 2
  fi

  perl -0pi -e "s/hsno_init = [0-9.]+_dbl_kind/hsno_init = ${thickness}_dbl_kind/" "${source_file}"
  echo "Set internal CICE hsno_init to ${thickness} m in ${source_file}"
}

case_prescribed_surface_flux() {
  local case_id="$1"

  case "${case_id}" in
    surface_warming)
      echo "${PRESCRIBED_WARMING_SURFACE_FLUX}"
      ;;
    surface_ablation)
      echo "${PRESCRIBED_ABLATION_SURFACE_FLUX}"
      ;;
    *)
      echo ""
      ;;
  esac
}

patch_prescribed_surface_flux() {
  local flux="$1"
  local source_file="${CICE_ROOT}/cicecore/cicedyn/general/ice_flux.F90"
  local replacement

  if [[ ! -f "${source_file}" ]]; then
    echo "CICE source file not found: ${source_file}" >&2
    return 2
  fi

  if ! grep -Eq "data fsurfn_d" "${source_file}"; then
    echo "Could not find fsurfn_d assignment in ${source_file}" >&2
    return 2
  fi

  replacement=$(printf "data fsurfn_d    / %s_dbl_kind, %s_dbl_kind, %s_dbl_kind, &\\n                          %s_dbl_kind, %s_dbl_kind, %s_dbl_kind /" \
                       "${flux}" "${flux}" "${flux}" "${flux}" "${flux}" "${flux}")
  SOURCE_FILE="${source_file}" REPLACEMENT="${replacement}" perl -0pi -e \
    's/data fsurfn_d\s+\/.*?\//$ENV{REPLACEMENT}/s' "${source_file}"
  echo "Set CICE prescribed fsurfn_d values to ${flux} W/m^2 in ${source_file}"
}

patch_balanced_bottom_flux_hook() {
  local enabled="$1"
  local source_file="${CICE_ROOT}/icepack/columnphysics/icepack_therm_vertical.F90"
  local flag=".false."

  if [[ "${enabled}" == "1" ]]; then
    flag=".true."
  fi

  if [[ ! -f "${source_file}" ]]; then
    echo "CICE source file not found: ${source_file}" >&2
    return 2
  fi

  if ! grep -q "climaseaice_balance_bottom_flux" "${source_file}"; then
    SOURCE_FILE="${source_file}" perl -0pi -e '
      s/character\(len=\*\),parameter :: subname='\''\(thermo_vertical\)'\''/logical (kind=log_kind), parameter :: climaseaice_balance_bottom_flux = .false.\n\n      character(len=*),parameter :: subname='\''(thermo_vertical)'\''/;
      s/fadvocn, saltvol, dfsalt ! advective heat flux to ocean/fadvocn, saltvol, dfsalt, \& ! advective heat flux to ocean\n         fbot_for_thickness      ! validation bottom flux for fixed-thickness cold case/;
      s/fadvocn = c0/fadvocn = c0\n      fbot_for_thickness = fbot/;
      s/(if \(icepack_warnings_aborted\(subname\)\) return\n\n         endif ! ktherm)/$1\n\n      if (climaseaice_balance_bottom_flux) fbot_for_thickness = fcondbotn/;
      s/smice,       massice,   \&\n                             smliq,       massliq,   \&\n                             fbot,        Tbot/smice,       massice,   \&\n                             smliq,       massliq,   \&\n                             fbot_for_thickness, Tbot/;
      s/fadvocn,   fbot/fadvocn,   fbot_for_thickness/;
    ' "${source_file}"
  fi

  if ! grep -q "climaseaice_balance_bottom_flux" "${source_file}"; then
    echo "Could not install balanced bottom-flux hook in ${source_file}" >&2
    return 2
  fi

  perl -0pi -e "s/climaseaice_balance_bottom_flux = \\.(?:true|false)\\./climaseaice_balance_bottom_flux = ${flag}/" "${source_file}"
  echo "Set CICE balanced bottom-flux hook to ${flag} in ${source_file}"
}

if [[ ! -x "${CICE_ROOT}/cice.setup" ]]; then
  echo "CICE setup script not found or not executable: ${CICE_ROOT}/cice.setup" >&2
  exit 2
fi

echo "CICE_ROOT=${CICE_ROOT}"
echo "CASE_ROOT=${CASE_ROOT}"
echo "MACHINE=${MACHINE}"
echo "ENVIRONMENT=${ENVIRONMENT}"
echo "PES=${PES}"
echo "RUN_CICE_CASES_EXECUTE=${EXECUTE}"
echo "RUN_CICE_CASES_BUILD=${BUILD}"
echo "RUN_CICE_CASES_RUN=${RUN}"
echo "CICE_INSTALL_LOCAL_TEMPLATES=${INSTALL_TEMPLATES}"
echo "CICE_CASE_FILTER=${CASE_FILTER:-<all>}"
echo "CICE_CONDUCT_FILTER=${CONDUCT_FILTER:-<all>}"
echo "CICE_CASE_SUFFIX=${CASE_SUFFIX:-<none>}"
echo "CICE_INTERNAL_SNOW_THICKNESS=${INTERNAL_SNOW_THICKNESS:-<unchanged>}"
echo "CICE_BALANCE_COLD_BOTTOM_FLUX=${BALANCE_COLD_BOTTOM_FLUX}"
echo "Case specs: ${REPO_ROOT}/validation/cice_bitz_lipscomb/cice_case_specs.toml"

mkdir -p "${CASE_ROOT}"

if [[ "${EXECUTE}" == "1" && "${INSTALL_TEMPLATES}" == "1" ]]; then
  install_local_templates
fi

if [[ "${EXECUTE}" == "1" && -n "${INTERNAL_SNOW_THICKNESS}" ]]; then
  patch_internal_snow_thickness "${INTERNAL_SNOW_THICKNESS}"
fi

for case_id in "${cases[@]}"; do
  if [[ -n "${CASE_FILTER}" && "${case_id}" != "${CASE_FILTER}" ]]; then
    continue
  fi

  for conduct in "${conducts[@]}"; do
    if [[ -n "${CONDUCT_FILTER}" && "${conduct}" != "${CONDUCT_FILTER}" ]]; then
      continue
    fi

    case_name="bl99_${case_id}_${conduct}${CASE_SUFFIX}"
    case_path="${CASE_ROOT}/${case_name}"
    case_runtime_set="$(runtime_set "${case_id}")"
    case_nml_set="$(case_set "${case_id}")"
    case_conduct_set="$(conduct_set "${conduct}")"
    case_surface_flux="$(case_prescribed_surface_flux "${case_id}")"
    setup_cmd=(
      "${CICE_ROOT}/cice.setup"
      --case "${case_path}"
      --mach "${MACHINE}"
      --env "${ENVIRONMENT}"
      --grid col
      --pes "${PES}"
      --set "ionetcdf,${case_runtime_set},histall,precision8,bl99_validation,${case_nml_set},${case_conduct_set}"
      --ignore-user-set
    )

    echo
    echo "==> ${case_name}"
    printf 'Setup command:'
    printf ' %q' "${setup_cmd[@]}"
    echo

    if [[ "${EXECUTE}" == "1" ]]; then
      if [[ "${BALANCE_COLD_BOTTOM_FLUX}" == "1" ]]; then
        if [[ "${case_id}" == "cold_conductive_relaxation" ]]; then
          patch_balanced_bottom_flux_hook 1
        else
          patch_balanced_bottom_flux_hook 0
        fi
      elif [[ -f "${CICE_ROOT}/icepack/columnphysics/icepack_therm_vertical.F90" ]] &&
           grep -q "climaseaice_balance_bottom_flux" "${CICE_ROOT}/icepack/columnphysics/icepack_therm_vertical.F90"; then
        patch_balanced_bottom_flux_hook 0
      fi

      if [[ -n "${case_surface_flux}" ]]; then
        patch_prescribed_surface_flux "${case_surface_flux}"
      fi

      (cd "${CICE_ROOT}" && "${setup_cmd[@]}")
      echo "Created ${case_path}"

      if [[ "${BUILD}" == "1" ]]; then
        (cd "${case_path}" && ./cice.build)
      fi

      if [[ "${RUN}" == "1" ]]; then
        (cd "${case_path}" && ./cice.run)
      fi
    else
      echo "Dry run only. Set RUN_CICE_CASES_EXECUTE=1 to create cases."
    fi
  done
done

if [[ "${EXECUTE}" == "1" && "${BALANCE_COLD_BOTTOM_FLUX}" == "1" ]]; then
  patch_balanced_bottom_flux_hook 0
fi
