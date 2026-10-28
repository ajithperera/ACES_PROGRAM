# Sourced by the cmake/*.sh build and test steps.
#
# aces_load_modules "<module list>": `module purge; module load <list>` when
# a list is given and environment modules exist (HiPerGator's default set
# is ACES_MODULES in ../CMakeLists.txt). With an empty list
# (cmake -DACES_MODULES="") the current environment is used as is: the
# Intel compilers (ifx/icx), MPI wrappers, MKL (MKLROOT) and python3 must
# then already be on PATH.
aces_load_modules() {
  local mods="$1" opts="$-"
  [ -z "$mods" ] && return 0
  set +u
  if ! type module >/dev/null 2>&1 && [ -f /etc/profile.d/modules.sh ]; then
    source /etc/profile.d/modules.sh
  fi
  if type module >/dev/null 2>&1; then
    module purge
    # shellcheck disable=SC2086
    module load $mods
  else
    echo "ACES_MODULES='$mods' but there is no 'module' command;" \
         "using the current environment" >&2
  fi
  case "$opts" in *u*) set -u ;; esac
}
