# EndoMT 가상 시간 기반 군집 색상 시각화 발표 자료

## 코드

```r
# 1. 시각화용 데이터 프레임(plot_df) 내 군집(Cluster) 메타데이터 매핑
# (각 세포의 고유 ID를 기준으로 소속 클러스터 정보를 데이터프레임에 병합합니다.)
plot_df$Cluster <- cell_clusters[rownames(plot_df)]

# 2. 가상 시간(Pseudotime)에 따른 주요 지표별 변화 추이 시각화
# (산점도 위에 LOESS 회귀선을 중첩하여 전반적인 발현 경향성을 분석합니다.)

# [A] Endothelial Score: 정상 내피세포 특성 지표의 감소 추세 분석
p_score_endo <- ggplot(plot_df, aes(x = Pseudotime, y = Endo_Score)) +
  geom_point(aes(color = Cluster), alpha = 0.5, size = 1) +
  geom_smooth(color = "black", method = "loess", size = 1.5) +
  theme_minimal() + ggtitle("Endothelial Score ( ↓ )")

# [B] Mesenchymal Score: 중간엽 특성 지표의 증가 추세 분석
p_score_mes <- ggplot(plot_df, aes(x = Pseudotime, y = Mes_Score)) +
  geom_point(aes(color = Cluster), alpha = 0.5, size = 1) +
  geom_smooth(color = "black", method = "loess", size = 1.5) +
  theme_minimal() + ggtitle("Mesenchymal Score ( ↑ )")

# [C] VWF Expression: 핵심 내피세포 마커 유전자의 발현 감소 분석
p_vwf <- ggplot(plot_df, aes(x = Pseudotime, y = VWF)) +
  geom_point(aes(color = Cluster), alpha = 0.5, size = 1) +
  geom_smooth(color = "black", method = "loess", size = 1.5) +
  theme_minimal() + ggtitle("VWF Expression ( ↓ )")

# [D] FN1 Expression: 핵심 EndoMT 마커 유전자의 발현 증가 분석
p_fn1 <- ggplot(plot_df, aes(x = Pseudotime, y = FN1)) +
  geom_point(aes(color = Cluster), alpha = 0.5, size = 1) +
  geom_smooth(color = "black", method = "loess", size = 1.5) +
  theme_minimal() + ggtitle("FN1 Expression ( ↑ )")

# 3. patchwork 패키지를 활용한 2x2 다중 패널(Multi-panel) 레이아웃 출력
(p_score_endo | p_score_mes) / (p_vwf | p_fn1)
```

## 발표 대본

이 코드는 앞에서 생성한 `plot_df`를 바탕으로, 가상 시간의 흐름에 따라 주요 지표들이 어떻게 변화하는지를 군집 정보까지 함께 반영하여 시각화하기 위한 코드입니다. 목적은 pseudotime 축을 따라 정상 내피세포 특성을 나타내는 점수와 유전자는 감소하는지, 반대로 EndoMT와 관련된 중간엽 특성 점수와 유전자는 증가하는지를 확인하는 동시에, 이러한 변화가 어떤 클러스터에 속한 세포들에서 주로 나타나는지를 함께 보여 주는 것입니다. 먼저 `plot_df$Cluster <- cell_clusters[rownames(plot_df)]`를 통해 각 세포의 고유 ID를 기준으로 군집 정보를 `plot_df`에 추가합니다. 이 단계는 단순히 발현값과 pseudotime만 보는 것이 아니라, 각 점이 어느 군집에 속하는 세포인지를 색으로 구분해서 해석할 수 있도록 만드는 과정입니다.

그다음 네 개의 그래프를 각각 생성하여 pseudotime에 따른 주요 지표의 변화를 시각화합니다. 첫 번째 그래프인 `p_score_endo`는 `Endo_Score`를 y축에 두고 `Pseudotime`을 x축에 두어, 정상 내피세포 특성이 가상 시간의 흐름에 따라 어떻게 변하는지를 보여 줍니다. 여기서 `geom_point(aes(color = Cluster), alpha = 0.5, size = 1)`는 각 세포를 점으로 표시하면서 소속 군집에 따라 색을 다르게 지정하는 부분입니다. `alpha = 0.5`는 점의 투명도를 조절하여 겹침이 심한 구간도 비교적 잘 보이도록 하고, `size = 1`은 점 크기를 지정합니다. 이어서 `geom_smooth(color = "black", method = "loess", size = 1.5)`는 전체 세포 분포를 바탕으로 LOESS 회귀 곡선을 검은색으로 그려서, 점들의 개별 변동과 별도로 전반적인 추세를 한눈에 볼 수 있게 합니다. 즉, 이 그래프는 군집별 분포와 전체 감소 경향을 동시에 보여 주는 구조입니다.

두 번째 그래프인 `p_score_mes`는 같은 방식으로 `Mes_Score`의 변화를 시각화합니다. 여기서는 중간엽 특성을 대표하는 점수가 pseudotime이 증가할수록 전반적으로 상승하는지를 확인하는 것이 핵심입니다. 세 번째 그래프인 `p_vwf`는 대표적인 내피세포 마커 유전자인 `VWF`의 발현 변화를 보여 주고, 네 번째 그래프인 `p_fn1`은 대표적인 EndoMT 또는 중간엽 관련 마커인 `FN1`의 발현 변화를 보여 줍니다. 따라서 이 네 개의 그래프를 함께 보면, 가상 시간의 증가에 따라 내피 특성은 감소하고 중간엽 특성은 증가하는 전환 패턴이 점수 수준과 개별 유전자 수준에서 동시에 관찰되는지를 확인할 수 있습니다. 또한 각 점이 군집별 색으로 표시되기 때문에, 특정 구간의 pseudotime에 어떤 클러스터가 주로 분포하는지도 함께 해석할 수 있습니다.

마지막으로 `(p_score_endo | p_score_mes) / (p_vwf | p_fn1)`를 사용하여 네 개의 그래프를 2행 2열 형태로 한 화면에 배치합니다. 위쪽에는 점수 기반 지표인 Endothelial Score와 Mesenchymal Score를, 아래쪽에는 대표 유전자 기반 지표인 VWF와 FN1을 배치함으로써, 세포 상태 전환을 점수와 개별 마커 두 수준에서 동시에 비교할 수 있도록 구성한 것입니다. 정리하면, 이 코드는 pseudotime 축을 따라 Endothelial Score, Mesenchymal Score, VWF, FN1의 변화를 군집 정보와 함께 시각화하여, EndoMT 과정이 단순한 군집 차이가 아니라 연속적인 상태 변화로 나타난다는 점을 보다 설득력 있게 제시하기 위한 코드입니다. 핵심 의의는 전체 추세선과 개별 세포 분포, 그리고 군집 정보를 동시에 제시함으로써, 세포 상태 전환의 방향성과 그 전환에 참여하는 군집 구조를 함께 설명할 수 있다는 점에 있습니다.
