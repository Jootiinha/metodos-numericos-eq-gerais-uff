# Métodos Numéricos - Equação de Helmholtz 2D

Este projeto implementa uma solução numérica da equação de Helmholtz 2D usando diferenças finitas centradas de 2ª ordem. O método resolve o problema em um domínio quadrado (0,1)×(0,1) com condições de contorno de Dirichlet impostas pela solução exata (combinação de 3 ondas planas).

## 📋 Sobre o Projeto

O projeto foi desenvolvido como parte do trabalho de Métodos Numéricos para Equações Gerais da UFF. Implementa:

- **Método**: Diferenças finitas centradas de 2ª ordem (stencil de 5 pontos)
- **Domínio**: (0,1) × (0,1) com malha uniforme
- **Condições de contorno**: Dirichlet (impostas pela solução exata)
- **Solução exata**: Combinação de 3 ondas planas com diferentes ângulos
- **6 grupos de configurações**: Cada grupo com diferentes combinações de ângulos

## 🚀 Instalação

### Pré-requisitos

- Python 3.11 ou superior
- pip ou poetry

### Passos

1. **Clone o repositório:**
   ```bash
   git clone <url-do-repositorio>
   cd metodos-numericos-eq-gerais-uff
   ```

2. **Instalar dependências:**
   ```bash
   # Usando pip
   pip install -r requirements.txt
   
   # Ou usando poetry
   poetry install
   ```

3. **Criar arquivo de configuração:**
   ```bash
   cp config_exemplo.yaml config.yaml
   # Edite config.yaml conforme necessário
   ```

## Uso

### Comando Recomendado
```bash
./run.sh  # Executa com configuração padrão (config.yaml)
```

### Otimizações de Desempenho
```bash
# Execução paralela (4 processos - recomendado para sistemas multi-core)
./run.sh --parallel 4

# Desabilitar cache (apenas para debug)
./run.sh --no_cache

# Combinar otimizações
./run.sh --parallel 4 --no_plots
```

### Outras Opções
```bash
# Usar configuração personalizada
./run.sh --config meu_config.yaml

# Sobrescrever parâmetros específicos
./run.sh --groups 1 2 --Ns 64 128 --ks 1 20 --no_plots

# Ver ajuda completa
./run.sh --help
```

> 📖 **Veja [OTIMIZACOES.md](OTIMIZACOES.md) para detalhes sobre otimizações de desempenho**

## Configuração

Todos os parâmetros podem ser configurados no arquivo `config.yaml`:

### Parâmetros de Simulação
- `groups`: Lista de grupos a executar (1-6)
- `Ns`: Tamanhos de malha (N subintervalos)
- `ks`: Valores do parâmetro k (número de onda)
- `solver`: Método de solução linear (`auto`, `splu`, `spsolve`, `gmres_ilu`)
- `plot_what`: Tipo de plot (`real` ou `abs`)
- `no_plots`: Desabilitar geração de figuras
- `group_thetas`: Definição dos ângulos para cada grupo

### Otimizações de Desempenho
- `parallel`: Número de processos paralelos (1=sequencial, recomendado: número de CPUs)
- `no_cache`: Desabilitar cache de matrizes (false=habilitado, recomendado: false)

## 📁 Estrutura do Projeto

```
.
├── main.py                 # Script principal
├── config.yaml            # Arquivo de configuração (criar a partir de config_exemplo.yaml)
├── config_exemplo.yaml    # Exemplo de configuração
├── run.sh                 # Script de execução conveniente
├── requirements.txt       # Dependências Python
├── pyproject.toml         # Configuração Poetry
│
├── outputs/               # Arquivos de saída (métricas CSV)
│   └── metrics_*.csv
│
├── figures/               # Gráficos gerados (PDF)
│   └── group_[1-6]/
│
├── slide/                 # Código LaTeX da apresentação
├── slide_overleaf/        # Versão para Overleaf
│
└── fontes_de_info/        # Material de referência (PDFs)
```

## 📊 Arquivos de Saída

- `outputs/metrics_all_groups.csv`: Métricas consolidadas de todos os grupos
- `outputs/metrics_group[1-6].csv`: Métricas individuais por grupo
- `figures/group_[1-6]/`: Gráficos 3D e cortes em PDF (se `no_plots: false`)

## ⚙️ Configuração

Veja `config_exemplo.yaml` para um exemplo completo de configuração. O arquivo `config.yaml` (não versionado) deve ser criado localmente.

### Parâmetros Principais

- **groups**: Lista de grupos a executar (1-6)
- **Ns**: Tamanhos de malha (ex: [64, 128, 192, 256])
- **ks**: Valores do parâmetro k/número de onda (ex: [1, 20, 40, 100])
- **solver**: Método de solução (`auto`, `splu`, `spsolve`, `gmres_ilu`)
- **plot_what**: Tipo de plot (`real` ou `abs`)
- **no_plots**: Desabilitar geração de figuras (útil para execuções rápidas)
- **parallel**: Número de processos paralelos (recomendado: número de CPUs)

## 📚 Documentação Adicional

- [OTIMIZACOES.md](OTIMIZACOES.md) - Detalhes sobre otimizações de desempenho
- [MELHORIAS_IMPLEMENTADAS.md](MELHORIAS_IMPLEMENTADAS.md) - Histórico de melhorias

## 🔬 Métodos Implementados

O projeto utiliza diferentes solvers lineares para resolver o sistema esparso:

- **auto**: Seleção automática baseada no tamanho do problema
- **splu**: Decomposição LU esparsa (eficiente para problemas médios)
- **spsolve**: Solver direto (para problemas pequenos)
- **gmres_ilu**: GMRES com pré-condicionador ILU (para problemas grandes)

## 📝 Notas

- Os arquivos de saída (`outputs/` e `figures/`) são ignorados pelo git (veja `.gitignore`)
- Crie seu próprio `config.yaml` a partir de `config_exemplo.yaml`
- Para execuções rápidas sem gráficos, use `--no_plots` ou `no_plots: true`

## 👤 Autor

João Carlos Romero Monteiro

## 📄 Licença

Este projeto foi desenvolvido para fins acadêmicos.
