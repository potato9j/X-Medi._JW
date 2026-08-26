# 초기 데이터 로드 및 QC-UMAP 분석 발표 자료

## 코드

```r
# 1. 필수 분석 라이브러리 로드
library(Seurat)

# 2. 데이터 경로 설정 및 샘플 식별자(Identity) 정의
# 분석 대상 데이터가 저장된 폴더 경로를 지정하고 각 샘플의 이름을 설정합니다.
main_dir <- choose.dir(caption = "Select the GSE162631_raw_counts_matrix directory")
main_dir <- gsub("\\", "/", main_dir)

sample_names <- c("R1_N", "R1_T", "R2_N", "R2_T", "R3_N", "R3_T", "R4_N", "R4_T")
paths <- file.path(main_dir, sample_names)
names(paths) <- sample_names

# 3. 10X Genomics 데이터 로드 및 Seurat 객체 생성
# 원본 데이터를 읽어 들여 'GBM_Project'라는 이름의 통합 분석 객체를 생성합니다.
raw_data <- Read10X(data.dir = paths)
merged_seurat <- CreateSeuratObject(counts = raw_data, project = "GBM_Project")

# 4. 메타데이터 정제: 바코드로부터 정확한 샘플 ID 추출
# 각 세포의 바코드 정보에서 샘플 출처를 식별할 수 있는 이름을 추출하여 저장합니다.
merged_seurat$orig.ident <- sub("^([^_]+_[^_]+)_.*", "\1", colnames(merged_seurat))
Idents(merged_seurat) <- "orig.ident"

# =====================================================================
# 데이터 품질 관리(Quality Control, QC) 및 필터링
# =====================================================================

# 5. 미토콘드리아 유전자 비율 계산
# 세포 사멸(Dead cells) 여부를 판단하기 위해 미토콘드리아 유전자 발현 비율을 산출합니다.
merged_seurat[["percent.mt"]] <- PercentageFeatureSet(merged_seurat, pattern = "^MT-")

# 6. 저품질 세포 제거 (Filtering)
# 유전자 검출 수가 너무 적거나(500개 미만) 너무 많은 경우(6000개 초과), 
# 그리고 미토콘드리아 비율이 높은(15% 초과) 저품질 세포를 분석에서 제외합니다.
merged_seurat <- subset(merged_seurat, subset = nFeature_RNA > 500 & nFeature_RNA < 6000 & percent.mt < 15)

# =====================================================================
# 데이터 정규화 및 선형 차원 축소
# =====================================================================

# 7. 데이터 정규화 및 고변동 유전자 식별
# 샘플 간 편차를 줄이기 위해 데이터를 정규화하고, 분석에 중요한 상위 2000개 유전자를 선정합니다.
merged_seurat <- NormalizeData(merged_seurat)
merged_seurat <- FindVariableFeatures(merged_seurat, selection.method = "vst", nfeatures = 2000)

# 8. 데이터 스케일링 및 주성분 분석(PCA) 실행
# 유전자 간 발현 단위를 맞추고, 데이터의 주요 특징을 추출하기 위한 선형 차원 축소(PCA)를 수행합니다.
merged_seurat <- ScaleData(merged_seurat)
merged_seurat <- RunPCA(merged_seurat, verbose = FALSE)

# =====================================================================
# 비선형 차원 축소 및 시각화 (UMAP)
# =====================================================================

# 9. 상위 30개 주성분을 활용한 UMAP 실행
# 고차원 데이터를 2차원 평면에 효과적으로 시각화하기 위해 UMAP 알고리즘을 적용합니다.
merged_seurat <- RunUMAP(merged_seurat, dims = 1:30)

# 10. 샘플 기원별 UMAP 시각화
# 각 세포가 어떤 샘플에서 왔는지 그룹화하여 최종 UMAP 지도를 출력합니다.
DimPlot(merged_seurat, reduction = "umap")
```

## 발표 대본

이 코드는 Glioblastoma 단일세포 RNA-seq 원시 데이터를 처음 불러오는 단계부터 품질 관리, 정규화, 차원 축소, 그리고 UMAP 시각화까지 수행하는 기본 전처리 파이프라인입니다. 목적은 여러 샘플에서 얻은 10X Genomics 원시 카운트 데이터를 하나의 Seurat 객체로 통합한 뒤, 품질이 낮은 세포를 제거하고, 분석 가능한 형태로 정리하여 각 샘플의 전반적인 분포를 UMAP 상에서 확인하는 것입니다. 먼저 `library(Seurat)`를 불러옵니다. Seurat는 단일세포 RNA-seq 데이터 분석에서 가장 기본이 되는 패키지로, 데이터 로드, 메타데이터 관리, 품질 관리, 정규화, 차원 축소, 시각화까지 전체 분석 흐름을 담당합니다.

그다음 `choose.dir()`를 사용하여 분석 대상 데이터가 들어 있는 상위 폴더를 선택하고, `gsub("\\", "/", main_dir)`를 통해 경로 구분자를 정리합니다. 이후 `sample_names <- c("R1_N", "R1_T", "R2_N", "R2_T", "R3_N", "R3_T", "R4_N", "R4_T")`로 총 8개 샘플의 이름을 지정하고, `file.path(main_dir, sample_names)`를 이용해 각 샘플 폴더의 전체 경로를 생성합니다. `names(paths) <- sample_names`는 각 경로에 샘플 이름을 연결하는 단계입니다. 이 부분의 의미는 여러 샘플의 원시 데이터를 한 번에 읽어 오기 위한 입력 구조를 만드는 것입니다. 이후 `Read10X(data.dir = paths)`를 통해 10X Genomics 형식의 raw count matrix를 불러오고, `CreateSeuratObject(counts = raw_data, project = "GBM_Project")`를 이용해 `merged_seurat`라는 Seurat 객체를 생성합니다. 즉, 원시 카운트 데이터를 단일세포 분석이 가능한 객체 형태로 변환하는 단계입니다.

이후 메타데이터를 정리합니다. `merged_seurat$orig.ident <- sub("^([^_]+_[^_]+)_.*", "\1", colnames(merged_seurat))`는 각 세포의 barcode 이름에서 샘플 출처를 나타내는 부분만 정규표현식으로 추출하여 `orig.ident`에 저장하는 코드입니다. 예를 들어 barcode 이름에 포함된 `R1_N`이나 `R2_T` 같은 샘플 정보를 분리해 내는 과정입니다. 이어서 `Idents(merged_seurat) <- "orig.ident"`를 설정하면, Seurat 객체에서 기본 identity가 샘플 출처 기준으로 바뀌게 됩니다. 이렇게 하면 이후 시각화나 비교 분석에서 세포를 샘플 단위로 쉽게 구분할 수 있습니다.

그다음 품질 관리 단계로 넘어갑니다. `PercentageFeatureSet(merged_seurat, pattern = "^MT-")`는 미토콘드리아 유전자의 발현 비율을 계산하여 `percent.mt`라는 메타데이터 열에 저장합니다. 미토콘드리아 유전자 비율은 손상되었거나 죽어가는 세포를 판별하는 대표적인 지표입니다. 이후 `subset(merged_seurat, subset = nFeature_RNA > 500 & nFeature_RNA < 6000 & percent.mt < 15)`를 사용하여 저품질 세포를 제거합니다. 여기서 `nFeature_RNA > 500`은 검출 유전자 수가 너무 적은 세포를 제거하는 조건이고, `nFeature_RNA < 6000`은 유전자 수가 과도하게 많은 세포를 제거하여 doublet 가능성을 줄이기 위한 조건입니다. 또 `percent.mt < 15`는 미토콘드리아 비율이 높은 세포를 제외하기 위한 기준입니다. 즉, 이 단계의 목적은 분석에 적합한 고품질 세포만 남기는 것입니다.

이후에는 데이터를 분석 가능한 형태로 정리합니다. `NormalizeData(merged_seurat)`는 세포마다 서로 다른 총 발현량 차이를 보정하여 비교 가능한 상태로 만들고, `FindVariableFeatures(merged_seurat, selection.method = "vst", nfeatures = 2000)`는 세포 간 차이를 잘 설명하는 상위 2000개의 고변동 유전자를 선택합니다. 그다음 `ScaleData(merged_seurat)`는 유전자 발현값의 스케일을 맞추어 특정 유전자가 과도하게 분석을 지배하지 않도록 조정하고, `RunPCA(merged_seurat, verbose = FALSE)`는 고차원 발현 데이터를 저차원으로 요약하는 주성분 분석을 수행합니다. PCA는 이후 UMAP 계산의 입력이 되는 핵심 단계입니다.

마지막으로 `RunUMAP(merged_seurat, dims = 1:30)`을 사용하여 상위 30개 주성분을 기반으로 UMAP을 계산합니다. UMAP은 고차원 데이터를 2차원 공간으로 축소하여 세포 간 유사성을 시각적으로 보여 주는 방법입니다. 이어서 `DimPlot(merged_seurat, reduction = "umap")`을 실행하면 최종 UMAP이 출력됩니다. 이 그래프를 통해 각 세포가 전체 공간에서 어떻게 분포하는지, 그리고 샘플별로 어떤 패턴을 보이는지를 확인할 수 있습니다. 정리하면, 이 코드는 여러 개의 10X Genomics 원시 데이터를 불러와 Seurat 객체를 생성하고, 샘플 메타데이터를 정리한 뒤, 품질 관리와 정규화, PCA, UMAP까지 수행하여 단일세포 데이터의 기본 구조를 시각화하는 초기 전처리 코드입니다. 핵심 의의는 이후의 통합 분석이나 세부 군집 분석에 들어가기 전에, 데이터 품질을 정제하고 전체 샘플 분포를 확인할 수 있는 기초 분석 기반을 마련한다는 점에 있습니다.
