# mhc-fold

 Arquivos do Marvin

└── marvin
    ├── afdb_up                                         # Alinhamento com AlphaFold/UniProt (Versão do target definitiva)
    │
    ├── afdb_up50                                       # Alinhamento com AlphaFold/UniProt_50 (Versão do target utilizada apenas para testes)
    │
    └── databases                                       # Databases utilizadas durante os alinhamentos


# Arquivos locais

├── query                                               # Estruturas de referência para o alinhamento com Foldseek
│   └── TCR3d
│       └── MHCI_groove_aligned
│           ├── db_filtrado                             # Dataset filtrado, inclui apenas estruturas de MHC HLA-A
│           └── db_samuel                               # Dataset completo, produzido por Samuel
│
│
└── version_02
    ├── alphafold                                       # Output de alinhamento do Foldseek p/ target AlphaFold
    │
    ├── pdb                                             # Output de alinhamento do Foldseek p/ target PDB
    │
    ├── analysis
    │   ├── binders                                     # Scripts de análises dos binders: anotação funcional, PCA e matriz de similaridade (tmalign)
    │   │   ├── functional_annotation
    │   │   ├── pca
    │   │   │   └── plots
    │   │   └── tmalign
    │   │
    │   ├── functional_annotation                       # Anotação funcional do MHC, p/ datasets AFDB e PDB
    │   │   ├── afdb
    │   │   ├── pdb                                     # Pastas "afdb" e "pdb" envolvem a coleta dos dados funcionais "brutos"
    │   │   │
    │   │   ├── filtered                                # Engloba os filtros globais: espécie humana, lista de estruturas em uniprot_to_remove, entre outros
    │   │   │   └── uniprot_to_remove
    │   │   │
    │   │   └── new_filters                             # Exclusivamente os filtros feitos para dividir AFDB em dois núcleos de análise distintos
    │   │
    │   ├──tmalign                                      # Várias versões da matriz de similaridade com diferentes combinações dos datasets
    │   │   ├── filter_length                           # Esses arquivos são utilizados no cálculo da PCA
    │   │   ├── filter_reviewed
    │   │   ├── filtered_afdb_first_version
    │   │   ├── filtered_pdb
    │   │   ├── pdb_vs_afdb
    │   │   └── reviewed+length
    │   │
    │   ├── pca                                         # Várias versões da PCA com diferentes combinações dos datasets
    │   │   ├── afdb_filter-length                      # As versões que entraram no relatório são: pdb, reviewed+length (AFDB com dois filtros juntos) e pdb_vs_afdb (combinação dos datasets)
    │   │   ├── afdb_filter-reviewed
    │   │   ├── afdb_first-version
    │   │   ├── pdb
    │   │   ├── pdb_vs_afdb
    │   │   └── reviewed+length
    │   │
    │   ├── helices                                     # Scripts, resultados e plots das análises de hélices
    │   │   ├── helix
    │   │   └── interface
    │   │       └── heatmap
    │   │
    │   ├── pykvfinder                                  # Análises de cavidade com pyKVFinder
    │   │   ├── afdb
    │   │   │   ├── afdb_aligned                        # Estruturas alinhadas com relação ao PDB de referência
    │   │   │   │   └── matrices
    │   │   │   └── afdb_cavities                       # Cavidades encontradas, uma pasta para cada estrutura. Contém as cavidades totais e filtradas
    │   │   │       └── AF-{uniprot_id}_mhc_aligned
    │   │   └── pdb
    │   │       ├── pdb_aligned
    │   │       │   └── matrices
    │   │       └── pdb_cavities
    │   │           └── {pdb_id}_MHC_groove_aligned
    │   │
    │   ├── resolution                                  # Scripts e plots da análise de resolução
    │   │
    │   └── plddt                                       # Scripts e plots das análises de pLDDT
    │
    └── filter                                          # Todas as estapas de tratamento estrutural
        ├── step0
        ├── step1                                       # Cada pasta de step{i} possui:
        ├── step2                                         # 1. O script que realiza o tratamento daquele step
        ├── step3                                         # 2. Os resultados divididos em duas pastas: pdb e afdb (quando aplicável)
        ├── step4                                         # 3. Em alguns casos, pastas adicionais com scripts de testes ou avaliações criados na programação do step
        ├── step5
        └── step6

# Criado por Giovanna de Lima em 04/04/2026
