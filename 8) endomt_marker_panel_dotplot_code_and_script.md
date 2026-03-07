# EndoMT 마커 패널 DotPlot 발표 자료

## PPT에 게시할 코드

```r
library(Seurat)
library(ggplot2)

# 1. 정확한 유전자 발현량 분석을 위해 기본 데이터(Assay)를 원본 RNA로 설정
DefaultAssay(endo_cells) <- "RNA"

# 2. EndoMT(내피-중간엽 이행) 현상 확인을 위한 핵심 마커 유전자 패널 구성
# - 앞쪽 3개: 정상 혈관 내피세포(Endothelial cell) 유지 마커
# - 뒤쪽 4개: 중간엽(Mesenchymal) 세포 및 형질 전환 마커
endoMT_markers <- c("PECAM1", "VWF", "CDH5", "FN1", "ACTA2", "TAGLN", "COL1A1")

# 3. DotPlot을 이용하여 세포 군집 간의 유전자 발현 패턴 변화(Switching) 시각화
DotPlot(endo_cells, features = endoMT_markers) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  ggtitle("Endothelial (↓) vs Mesenchymal (↑) Marker Panel")
```

## 발표 대본

이 코드는 혈관 내피세포 하위 군집 객체인 `endo_cells`를 대상으로, EndoMT 즉 내피-중간엽 이행 현상이 실제로 나타나는지를 핵심 마커 유전자 패널 수준에서 확인하기 위한 시각화 코드입니다. 목적은 정상 혈관 내피세포의 특성을 유지하는 마커와, 중간엽 세포 성격 또는 형질 전환과 관련된 마커를 한 그래프 안에서 동시에 비교하여, 특정 군집에서 내피 특성이 감소하고 중간엽 특성이 증가하는 발현 전환 양상이 나타나는지를 확인하는 것입니다. 먼저 `library(Seurat)`와 `library(ggplot2)`를 불러옵니다. `Seurat`는 단일세포 RNA-seq 데이터의 시각화와 군집 기반 발현 비교에 사용되는 핵심 패키지이고, `ggplot2`는 그래프 제목과 축 텍스트와 같은 시각적 요소를 조정하는 데 사용됩니다.

그다음 `DefaultAssay(endo_cells) <- "RNA"`를 실행하여 분석의 기준이 되는 assay를 원본 RNA 발현값으로 설정합니다. 이는 하위 군집 수준에서 특정 마커 유전자의 실제 발현량을 해석할 때, 통합본이나 다른 보정 결과가 아니라 원본 RNA 발현을 기준으로 비교하겠다는 의미입니다. 이후 `endoMT_markers <- c("PECAM1", "VWF", "CDH5", "FN1", "ACTA2", "TAGLN", "COL1A1")`를 통해 EndoMT 해석에 필요한 핵심 유전자 패널을 정의합니다. 여기서 `PECAM1`, `VWF`, `CDH5`는 정상 혈관 내피세포의 정체성을 유지하는 대표적인 마커이고, `FN1`, `ACTA2`, `TAGLN`, `COL1A1`은 중간엽 세포 특성 또는 형질 전환 상태와 관련된 마커입니다. 따라서 이 유전자 패널은 내피세포 상태와 중간엽 상태를 한 번에 비교하기 위한 구조로 설계되어 있습니다.

마지막으로 `DotPlot(endo_cells, features = endoMT_markers)`를 사용하여 각 군집별 유전자 발현 패턴을 점 그래프 형태로 시각화합니다. DotPlot에서 각 점의 크기는 해당 군집에서 그 유전자를 발현하는 세포 비율을 나타내고, 점의 색은 평균 발현량의 상대적 수준을 나타냅니다. 따라서 점이 크고 색이 진하면 해당 군집에서 많은 세포가 그 유전자를 높게 발현하고 있다는 뜻입니다. 이후 `theme(axis.text.x = element_text(angle = 45, hjust = 1))`는 x축 유전자 이름이 겹치지 않도록 45도 기울여 표시하는 설정이고, `ggtitle("Endothelial (↓) vs Mesenchymal (↑) Marker Panel")`은 발표용으로 그래프 제목을 추가한 부분입니다. 정리하면, 이 코드는 혈관 내피세포 하위 군집에서 내피 마커와 중간엽 마커를 동시에 비교하는 DotPlot을 생성하여, 특정 군집에서 Endothelial 특성은 감소하고 Mesenchymal 특성은 증가하는 EndoMT형 발현 전환 양상이 나타나는지를 확인하기 위한 코드입니다. 핵심 의의는 개별 유전자 하나만 보는 것이 아니라, 기능적으로 연결된 마커 패널 전체를 함께 비교함으로써 EndoMT 해석의 신뢰도를 높일 수 있다는 점에 있습니다.
