import re
import numpy as np


def load_data():
    text = open('io/loadData.m').read()
    v_match = re.search(r'data_V\s*=\s*\[(.*?)\];', text, re.S)
    jd_match = re.search(r'data_JD\s*=\s*\[(.*?)\];', text, re.S)
    V = np.fromstring(v_match.group(1).replace('\n', ' '), sep=' ')
    JD = np.fromstring(jd_match.group(1).replace('\n', ' '), sep=' ')
    return V, JD


def errorFunction_loop(predicted, data_V, data_JD):
    threshold = 1e-12
    max_val = max(1e-12, np.max(np.abs(data_JD)))
    err = np.zeros_like(data_JD)
    for i in range(len(data_JD)):
        a = abs(data_JD[i])
        p = abs(predicted[i])
        if a < threshold or p < threshold:
            err[i] = (predicted[i] - data_JD[i]) / max_val
        else:
            err[i] = np.log10(p) - np.log10(a)
            if np.sign(predicted[i]) != np.sign(data_JD[i]) and p > threshold and a > threshold:
                err[i] *= 3
        if abs(data_V[i]) > 0.05:
            err[i] *= (1 + abs(err[i]))
    return err


def errorFunction_vec(predicted, data_V, data_JD):
    threshold = 1e-12
    max_val = max(1e-12, np.max(np.abs(data_JD)))
    a = np.abs(data_JD)
    p = np.abs(predicted)
    mask_small = (a < threshold) | (p < threshold)
    err = np.zeros_like(data_JD)
    err[mask_small] = (predicted[mask_small] - data_JD[mask_small]) / max_val
    mask_log = ~mask_small
    err[mask_log] = np.log10(p[mask_log]) - np.log10(a[mask_log])
    sign_mismatch = (np.sign(predicted) != np.sign(data_JD)) & (p > threshold) & (a > threshold)
    err[sign_mismatch] *= 3
    idx = np.abs(data_V) > 0.05
    err[idx] = err[idx] * (1 + np.abs(err[idx]))
    return err


def errorFunctionPartial_loop(predicted, data_V, data_JD, prev=None):
    threshold = 1e-12
    max_val = max(1e-12, np.max(np.abs(data_JD)))
    err = np.zeros_like(data_JD)
    for i in range(len(data_JD)):
        a = abs(data_JD[i])
        p = abs(predicted[i])
        if a < threshold or p < threshold:
            err[i] = (predicted[i] - data_JD[i]) / max_val
        else:
            err[i] = np.log10(p) - np.log10(a)
            if np.sign(predicted[i]) != np.sign(data_JD[i]):
                err[i] *= 4
        if abs(data_V[i]) > 0.05:
            if prev is not None:
                err[i] *= (1 + abs(prev[i]))
            else:
                err[i] *= (1 + abs(err[i]))
    return err


def errorFunctionPartial_vec(predicted, data_V, data_JD, prev=None):
    threshold = 1e-12
    max_val = max(1e-12, np.max(np.abs(data_JD)))
    a = np.abs(data_JD)
    p = np.abs(predicted)
    mask_small = (a < threshold) | (p < threshold)
    err = np.zeros_like(data_JD)
    err[mask_small] = (predicted[mask_small] - data_JD[mask_small]) / max_val
    mask_log = ~mask_small
    err[mask_log] = np.log10(p[mask_log]) - np.log10(a[mask_log])
    sign_mismatch = (np.sign(predicted) != np.sign(data_JD))
    err[sign_mismatch & mask_log] *= 4
    idx = np.abs(data_V) > 0.05
    if prev is not None:
        err[idx] = err[idx] * (1 + np.abs(prev[idx]))
    else:
        err[idx] = err[idx] * (1 + np.abs(err[idx]))
    return err


# Enhanced positive shares the same algorithm as Partial
errorFunctionEnhanced_loop = errorFunctionPartial_loop
errorFunctionEnhanced_vec = errorFunctionPartial_vec


def main():
    V, JD = load_data()
    predicted = JD * (1 + 0.1 * np.sin(V))
    prev = 0.5 * np.abs(np.sin(V))
    e_loop = errorFunction_loop(predicted, V, JD)
    e_vec = errorFunction_vec(predicted, V, JD)
    assert np.allclose(e_loop, e_vec)

    e_loop_p = errorFunctionPartial_loop(predicted, V, JD, prev)
    e_vec_p = errorFunctionPartial_vec(predicted, V, JD, prev)
    assert np.allclose(e_loop_p, e_vec_p)

    e_loop_e = errorFunctionEnhanced_loop(predicted, V, JD, prev)
    e_vec_e = errorFunctionEnhanced_vec(predicted, V, JD, prev)
    assert np.allclose(e_loop_e, e_vec_e)

if __name__ == '__main__':
    main()
