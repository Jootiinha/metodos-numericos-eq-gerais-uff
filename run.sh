#!/bin/bash
# Script para executar o programa Helmholtz 2D

# Ativar ambiente virtual
source .venv/bin/activate

# Executar o programa com todos os argumentos passados
python main.py "$@"

# ./run.sh --groups 1 2 --Ns 64 --ks 1 20 --no_plots --parallel 2
# ./run.sh --parallel 4  # Cache já está habilitado por padrão