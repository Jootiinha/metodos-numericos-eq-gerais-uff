# Apresentação Beamer - Solução Numérica da Equação de Helmholtz 2D

Esta pasta contém a apresentação em LaTeX/Beamer sobre o trabalho de solução numérica da equação de Helmholtz 2D.

## Estrutura Modular

A apresentação está organizada em arquivos separados para facilitar edição:

- `apresentacao.tex`: Arquivo principal (inclui todos os outros)
- `preambulo.tex`: Preâmbulo, pacotes e configurações
- `introducao.tex`: Slides da seção de introdução
- `metodologia.tex`: Slides da seção de metodologia
- `resultados.tex`: Slides da seção de resultados
- `conclusoes.tex`: Slides da seção de conclusões
- `final.tex`: Slide final
- `Makefile`: Script para compilar a apresentação
- `README.md`: Este arquivo

## Vantagens da Estrutura Modular

- **Fácil edição**: Cada seção em arquivo separado
- **Organização**: Fácil localizar e modificar conteúdo específico
- **Colaboração**: Múltiplas pessoas podem editar seções diferentes
- **Manutenção**: Mudanças isoladas não afetam outras seções

## Compilação

### Usando Make (Recomendado)

```bash
# Compilar a apresentação
make

# Compilar e abrir automaticamente
make view

# Limpar arquivos auxiliares
make clean

# Limpar tudo (incluindo PDF)
make cleanall
```

### Compilação Manual

```bash
# Compilar (2 passagens para sumário)
pdflatex apresentacao.tex
pdflatex apresentacao.tex

# Ou compilação rápida (1 passagem)
pdflatex apresentacao.tex
```

## Requisitos

- LaTeX com Beamer (geralmente incluído em distribuições como TeX Live ou MiKTeX)
- Pacotes necessários:
  - `beamer`
  - `babel` (português)
  - `amsmath`, `amssymb`
  - `graphicx`
  - `booktabs`
  - `siunitx`

## Estrutura da Apresentação

1. **Introdução**: Equação de Helmholtz e motivação
2. **Metodologia**: Diferenças finitas e solvers
3. **Resultados**: Erros, convergência e desempenho
4. **Conclusões**: Resultados principais e trabalhos futuros

## Personalização

Para personalizar a apresentação:

- Edite `apresentacao.tex` diretamente
- Modifique o tema em `\usetheme{Madrid}`
- Adicione figuras na pasta `slide/` e referencie com `\includegraphics`
- Ajuste cores com `\usecolortheme`

## Notas

- A apresentação usa aspect ratio 16:9 (`aspectratio=169`)
- Todas as figuras mencionadas devem estar na mesma pasta ou ajustar caminhos
- Para incluir figuras geradas pelo código Python, copie os PDFs de `figures/` para `slide/`
