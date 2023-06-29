# Matlab을 통해서 통신공학의 여러가지 AM modulation 방식을 공부합니다.
## DSB-SC
   - DSB-SC에서 수신자는 캐리어 신호의 정보를 알고있다고 가정하고 Demodulation을 했습니다.
   - 이상적인 필터를 사용해서 Modulation과 Demodulation을 진행했습니다.
## DSB-LC
  - 수신자는 캐리어 신호의 정보를 사전에 몰라도 수신 신호에 포함되어 있는 캐리어 정보를 활용하여 Demodulation 합니다.
  - envelope detector를 활용해 Demodulation 과정을 진행했습니다.
## SSB-SC
- 위의 두 방식과 다르게 single sideband만을 이용해 주파수 대역을 효율적으로 사용했습니다.
- Hilbert transform을 활용하여 upper sideband만을 사용해 송신과정을 나타냈습니다.
## VSB
- 이상적인 필터를 사용하는 것이 아닌 좀더 현실적인 필터를 사용해 Modulation 을 했습니다.
- 필터의 모양은 사다리꼴 형태로 만들어 선형성과 일정한 상수값을 유지하도록 했습니다.
