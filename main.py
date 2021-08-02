import numpy
from scipy.integrate import odeint
import matplotlib.pyplot as plt


def seir(seir, t, b, a, g, N):
    """
    seir: seir[0]: S, seir[1]: E, seir[2]: I, seir[3]: R
    t: 時間
    b: 感染率
    a: 発症率
    g: 回復率
    N: 全人口
    """
    return [
        -b * seir[0] * seir[2] / N,
        b * seir[0] * seir[2] / N - a * seir[1],
        a * seir[1] - g * seir[2],
        g * seir[2]
    ]


if __name__ == '__main__':
    """
    """
    R0 = 1.77 # 基本生産数
    infectious = 4.8  # 平均発症期間
    latency = 5.6 # 平均潜伏期間
    S = 13960000   # 無免疫者数
    E = R0 * 2195  # 潜伏感染者数
    I = 2195  # 発症者数
    R = 0  # 有免疫者数
    N = S + E + I + R  # 全人口
    isolation = 0.1  # 自粛(孤立やマスク着用)
    b = (1 - isolation) * R0 / infectious  # 感染率
    a = 1 / latency  # 1日あたりの発症率
    g = 1 / infectious  # 1日あたりの回復率
    all_days = 30
    t = numpy.arange(0, all_days, 1)
    args = (b, a, g, N)
    result = odeint(seir, (S, E, I, R), t, args)
    for days, data in enumerate(result):
        print(days, data[2])
    plt.plot(t, result)
    plt.legend(['S', 'E', 'I', 'R'])
    plt.xlabel('days')
    plt.ylabel('population')
    plt.show()
