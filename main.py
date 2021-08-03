import numpy
from scipy.integrate import odeint
import matplotlib.pyplot as plt
import click


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


@click.command()
@click.argument('total_days', default=30)
@click.argument('R0', default=1.77)  # 基本再生産数
@click.argument('average_infectious_period', default=4.8)  # 平均発症期間
@click.argument('average_incubation_period', default=5.6)  # 平均潜伏期間
@click.argument('susceptible', default=13960000)  # 無免疫者数
@click.argument('exposed', default=3885.15)  # 潜伏感染者数
@click.argument('infectious', default=2195)  # 発症者数
@click.argument('recovered', default=0)  # 有免疫者数
@click.argument('isolation', default=0.0)  # 自粛などの隔離をしている割合(0~1)
def run(r0, total_days, average_infectious_period, average_incubation_period, susceptible, exposed, infectious, recovered, isolation):
    S = susceptible
    E = exposed
    I = infectious
    R = recovered
    N = S + E + I + R  # 全人口
    b = (1 - isolation) * r0 / average_infectious_period  # 感染率
    a = 1 / average_incubation_period  # 1日あたりの発症率
    g = 1 / average_infectious_period  # 1日あたりの回復率
    t = numpy.arange(0, total_days, 1)
    args = (b, a, g, N)
    result = odeint(seir, (S, E, I, R), t, args)
    for days, data in enumerate(result):
        print(days, data[2])
    plt.plot(t, result)
    plt.legend(['Susceptible', 'Exposed', 'Infectious', 'Recovered'])
    plt.xlabel('days')
    plt.ylabel('population')
    plt.show()


if __name__ == '__main__':
    run()
