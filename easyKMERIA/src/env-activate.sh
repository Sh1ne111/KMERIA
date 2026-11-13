# 激活软件环境: 自动找 micromamba / mamba / conda (不依赖 .bashrc 的初始化, 哪个装了这个环境就用哪个)
_easy_rc_val() {  # 只读取 ~/.bashrc 等文件里 NAME=路径 的写法, 不执行它们
  local k="$1" f v
  for f in "$HOME/.bashrc" "$HOME/.bash_profile" "$HOME/.profile"; do
    [ -r "$f" ] || continue
    v="$(grep -E "^[[:space:]]*(export[[:space:]]+)?${k}=" "$f" 2>/dev/null | tail -n 1 | sed -E "s/^[^=]*=//; s/^[\"']//; s/[\"'].*\$//")" || v=""
    v="${v/#\$\{HOME\}/$HOME}"; v="${v/#\$HOME/$HOME}"; v="${v/#\~/$HOME}"
    case "$v" in
      /*) case "$v" in *[!A-Za-z0-9_./+-]*) ;; *) echo "$v"; return 0 ;; esac ;;
    esac
  done
  return 1
}
_easy_find() {  # 用法: _easy_find 程序名 [已知路径]; 依次看 PATH、已知路径、~/.bashrc 里写过的路径、常见安装位置
  local n="$1" hint="${2:-}" d p f
  p="$(type -P "$n" 2>/dev/null)" || p=""
  if [ -n "$p" ]; then echo "$p"; return 0; fi
  if [ -n "$hint" ] && [ -x "$hint" ]; then echo "$hint"; return 0; fi
  for f in "$HOME/.bashrc" "$HOME/.bash_profile" "$HOME/.profile"; do
    [ -r "$f" ] || continue
    while IFS= read -r p; do
      if [ -f "$p" ] && [ -x "$p" ]; then echo "$p"; return 0; fi
    done < <(grep -ohE "/[A-Za-z0-9_./+-]*/${n}" "$f" 2>/dev/null | awk '!s[$0]++')
  done
  for d in "$HOME/micromamba/bin" "$HOME/.local/bin" "$HOME/bin" "$HOME/miniforge3/bin" "$HOME/mambaforge/bin" "$HOME/miniconda3/bin" "$HOME/anaconda3/bin" "$HOME/mamba/bin" "$HOME/conda/bin" "$HOME"/*/bin /opt/conda/bin /opt/miniforge3/bin /opt/miniconda3/bin /opt/anaconda3/bin /usr/local/bin /share/apps/conda/bin; do
    if [ -f "$d/$n" ] && [ -x "$d/$n" ]; then echo "$d/$n"; return 0; fi
  done
  return 1
}
_easy_activate() {
  local env="$1" mm="" mba="" cda="" root="" r pick="" base found=""
  mm="$(_easy_find micromamba "${MAMBA_EXE:-}")" || mm=""
  mba="$(_easy_find mamba)" || mba=""
  cda="$(_easy_find conda "${CONDA_EXE:-}")" || cda=""
  [ -n "$mm" ] && found="$found micromamba=$mm"
  [ -n "$mba" ] && found="$found mamba=$mba"
  [ -n "$cda" ] && found="$found conda=$cda"
  if [ -z "$found" ]; then
    echo "[提示] 没有找到 micromamba / mamba / conda, 已跳过环境激活。请在服务器终端运行 type -a micromamba conda 看它在哪,"
    echo "       再把 export PATH=\"那个目录:\$PATH\" 写进表单里的“额外环境命令”。(当前 PATH=$PATH)"
    return 1
  fi
  # micromamba 的环境放在 MAMBA_ROOT_PREFIX 下; 没设置时试 .bashrc 里写的值和几个常见位置
  if [ -n "$mm" ]; then
    for r in "${MAMBA_ROOT_PREFIX:-}" "$(_easy_rc_val MAMBA_ROOT_PREFIX)" "$(dirname "$(dirname "$mm")")" "$HOME/micromamba" "$HOME/.local/share/mamba" "$HOME/mamba" "$HOME/.mamba"; do
      if [ -n "$r" ] && [ -d "$r/envs/$env" ]; then root="$r"; pick="micromamba"; break; fi
    done
  fi
  if [ -z "$pick" ] && [ -n "$mba" ] && "$mba" env list 2>/dev/null | grep -Fwq -- "$env"; then pick="mamba"; fi
  if [ -z "$pick" ] && [ -n "$cda" ] && "$cda" env list 2>/dev/null | grep -Fwq -- "$env"; then pick="conda"; fi
  if [ -z "$pick" ]; then
    if [ -n "$mm" ]; then pick="micromamba"; elif [ -n "$mba" ]; then pick="mamba"; else pick="conda"; fi
  fi
  case "$pick" in
    micromamba)
      export MAMBA_EXE="$mm"
      [ -n "$root" ] && export MAMBA_ROOT_PREFIX="$root"
      eval "$("$mm" shell hook --shell bash)"; micromamba activate "$env" ;;
    *)
      if [ "$pick" = mamba ]; then base="$("$mba" info --base 2>/dev/null)"; else base="$("$cda" info --base 2>/dev/null)"; fi
      if [ -f "$base/etc/profile.d/conda.sh" ]; then source "$base/etc/profile.d/conda.sh"; conda activate "$env"
      else eval "$("${mba:-$cda}" shell hook --shell bash)"; "$pick" activate "$env"; fi ;;
  esac
  if [ "${CONDA_DEFAULT_ENV:-}" = "$env" ] || [ "${MAMBA_DEFAULT_ENV:-}" = "$env" ]; then
    echo "[环境] 已用 $pick 激活 $env"
  else
    echo "[警告] 用 $pick 激活 $env 失败 (找到的程序:$found)。请确认环境名正确, 或在“额外环境命令”里写你自己的激活命令"; return 1
  fi
}
