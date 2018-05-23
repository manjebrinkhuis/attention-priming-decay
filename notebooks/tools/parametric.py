import numpy as np

# Parametric models
def exp_decay(x, w, t):
    """ """
    x = np.array(x)
    return w * np.exp(-x*t)


# Parametric models
def exp_decay_norm(x, t):
    """ """
    x = np.array(x)
    return (np.exp(-x*t)) * np.exp(-x*t)


def repeat(arr):
    count = 0
    out = [count]
    for i, val in enumerate(arr[1:]):
        count = count + 1 if val == arr[i-1] else 0
        out.append(count)

    return np.array(out)


def iter_model(xs, t=0, start=0.0):
    y = [start]
    for x in xs[:-1]:
        y.append(y[-1] + (((x)-(y[-1]))*t))

    return np.array(y)


def diff_model(xs, t=1, y0=.5, c=1):
    y = [y0]
    for x in xs[:-1]:
        y1 = y[-1]
        y2 = abs(x*c-((x*c-y1) * np.exp(-t)))
        y.append(y2)

    return np.array(y)


def decay_model(xs, t=1, y0=.5):
    y = [y0]
    for x in xs[:-1]:
        y1 = y[-1]
        y2 = y1 * np.exp(-t)
        y.append(y2)

    return np.array(y)


def martini_decay(x, w1, t1, w2, t2):
    x = np.array(x)
    return w1*np.exp(-(x/t1)) + w2*np.exp(-(x/t2))
