#!/usr/bin/env bash
# 第一次运行会自动建虚拟环境并安装依赖
cd "$(dirname "$0")"
[ -d .venv ] || python3 -m venv .venv
. .venv/bin/activate
pip install -q -r requirements.txt
python app.py --open "$@"
