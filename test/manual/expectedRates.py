"""
Quick python script to find the expected WS and RS decay rates

"""
import math
import matplotlib.pyplot as plt
import numpy as np


def fit_params(decay_params):
    """
    Find a, b, c from a dict of decay params
    """
    a = decay_params["r"] ** 2
    b = (
        decay_params["r"]
        * decay_params["width"]
        * (
            decay_params["y"] * decay_params["z_re"]
            + decay_params["x"] * decay_params["z_im"]
        )
    )
    c = (
        0.25
        * (decay_params["x"] ** 2 + decay_params["y"] ** 2)
        * decay_params["width"] ** 2
    )

    return a, b, c


def rs_rate(decay_params, time):
    return math.exp(-1 * decay_params["width"] * time)


def ws_rate(decay_params, a, b, c, time):
    return (a + b * time + c * time ** 2) * math.exp(-1 * decay_params["width"] * time)


def main():
    decay_params = {
        "x": 0.004,
        "y": 0.007,
        "r": 0.05,
        "z_im": -0.3,
        "z_re": 0.8,
        "width": 2500.0,
    }
    a, b, c = fit_params(decay_params)

    times = np.linspace(0, 0.002)
    rs_rates = [rs_rate(decay_params, time) for time in times]
    ws_rates = [ws_rate(decay_params, a, b, c, time) for time in times]

    plt.plot(times, rs_rates, label="RS")
    plt.legend()
    plt.show()

    plt.plot(times, ws_rates, label="WS")
    plt.legend()
    plt.show()


if __name__ == "__main__":
    main()
