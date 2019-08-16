# coding: utf-8
"""
# How to run??
ex.1) python NeuralNetwork.py --func "y=sin(x)" --lambda_func "lambda x:np.sin(x)" --xmin -3 --xmax 3
ex.2) python NeuralNetwork.py --func "y=|x|"    --lambda_func "lambda x:abs(x)" --N 30  --M 5
ex.3) python NeuralNetwork.py --func "y=x^2"    --lambda_func "lambda x:x*x" 
"""

import argparse
import numpy as np
from scipy import linalg as LA
import matplotlib.pyplot as plt

def forward(x, w1, w2):
    """
    関数の概要　：順伝播を計算する。
    @param x  ：入力
    @param w1 ：隠れ層の重み w1[i, j] はj番目の入力と隠れ層のi番目の素子の間の重み。
    @param w2 ：出力層の重み w2[i, j] は隠れ層のj番目素子と出力層のi番目の素子の間の重み。
    @return a1：a1[i]: 隠れ層iへの入力
    @return a2：a2[i]: 出力層iへの入力
    """
    a1 = w1.dot(np.append(x, 1))           # 隠れ層の入力を計算
    a2 = w2.dot(np.append(np.tanh(a1), 1)) # 出力層の入力を計算
    return (a1, a2)

def backprop(a1, a2, w1, w2, delta2):
    """
    関数の概要：誤差逆伝播を計算する。 h'(a) * Σw*δ
    @param a1, a2：各層への入力
    @param w1, w2：各層の重み(w2[:,0:M]で、バイアス項の重みを除く。)
    @param delta2：出力層の誤差
    @return 隠れ層の誤差
    """ # 逆伝播の場合は、重みのベクトルを転地する必要があることに注意！！
    return ((1- np.tanh(a1)**2)*(w2[:,0:M].T).dot(delta2)) # 隠れ層の誤差

# 偏微分係数の計算
def diffcoef(x, a1, a2, w1, w2, delta2):
    """
    関数の概要：偏微分係数を計算する
    @param x     ：入力
    @param a1,a2 ：隠れ層、出力層への入力
    @param w1,w2 ：隠れ層・出力層の重み
    @param delta2：出力結果の誤差
    """
    delta1 = backprop(a1, a2, w1, w2, delta2)
    diff1 = np.outer(delta1, np.append(x, 1))
    diff2 = np.outer(delta2, np.append(np.tanh(a1), 1))
    return (diff1, diff2)

def steepest_step(x, t, w1, w2):
    """
    関数の概要：最急降下法のステップ方向を計算する。
    @param x    ：入力
    @param t    ：正解
    @param w1,w2：隠れ層・出力層の重み
    @return diff：ステップ方向
    """
    diff1 = np.zeros((M, D+1))
    diff2 = np.zeros((K, M+1))
    for i in range(N):
        # i番目の学習データを順伝播させると、出力結果が返ってくる。
        # a1:各隠れ層(M個)への入力ベクトル。 a2:各出力層(K個)への入力ベクトル。
        a1, a2 = forward(x[i], w1, w2)
        # i番目の学習データに対する誤差関数の a1, a2 での偏微分係数を求める。
        d1, d2 = diffcoef(x[i], a1, a2, w1, w2, a2-t[i])
        diff1 -= ALPHA*d1 # 最急降下法なので、勾配の逆方向に進む。
        diff2 -= ALPHA*d2 # 最急降下法なので、勾配の逆方向に進む。
    return (diff1, diff2)

def steepest_descent_method(x, t):
    """
    関数の概要：最急降下法を用いて重みを最適化する。
    @param x：入力
    @param t：正解
    @return w1,w2：最適な重み
    @return i+1  ：反復回数
    """
    # 重み（行:次の層の数、列:今の層の数）の初期化
    w1 = np.random.uniform(-1, 1, (M, D+1))
    w2 = np.random.uniform(-1, 1, (K, M+1))
    for i in range(ITER_MAX):
        d1, d2 = steepest_step(x, t, w1, w2)
        w1 += d1
        w2 += d2
        # 共に誤差が十分小さくなったら終了する。
        if LA.norm(d1) < ITER_EPS and LA.norm(d2) < ITER_EPS:
            break
    return (w1, w2, i+1)

#=== フィッティング ===
def fit(func, lambda_func, xmin=-1, xmax=1):
    """
    関数の概要：与えられた関数をニューラルネットワークで学習する
    @param func       ：関数の式
    @param lambda_func：関数の機能
    @param xmin/xmax  ：定義域
    """
    #=== 正解の関数を図示 ===
    x = np.linspace(xmin, xmax, N)
    t = np.vectorize(lambda_func)(x)
    plt.xlim(xmin, xmax)
    plt.scatter(x, t, color="blue", label="learning data")
    plt.plot(x, t, color="blue", label="target function")

    #=== 重みの最適化(最急降下法) ===
    w1, w2, count = steepest_descent_method(x, t)
    y = np.vectorize(lambda x: forward(x, w1, w2)[1][0])(x)
    plt.plot(x, y, color="red", label="result")
    plt.title("{} (iteration={}".format(func, count))
    plt.legend()
    plt.savefig('{}.png'.format(func)) # --

    # === 中間層の結果(基底ベクトル)を表示する。
    plt.clf()
    for i in range(M):
        phi = np.vectorize(lambda x: np.tanh(forward(x, w1, w2)[0][i]))(x)
        plt.plot(x, phi, label="phi_{}".format(i+1))
    plt.scatter(x, t, color="blue", label="learning data")
    plt.title("{} (iteration={}".format(func, count))
    plt.legend()
    plt.savefig('{}-hidden-units.png'.format(func)) # --

def create_lambda_with_globals(s):
    return eval(s, globals())

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('--func', type=str, help='specify the function. (ex. $y=\sin(\pi x)$)' , required=True)
    parser.add_argument('--lambda_func', type=create_lambda_with_globals, help="Define the func as lambda function. (ex.'lambda x:np.sin(x)')", required=True)
    parser.add_argument('--xmin', type=float, default=-1, help="Define the minimum x in domain.")
    parser.add_argument('--xmax', type=float, default=1,  help="Define the maximum x in domain.")
    parser.add_argument('--N', type=int, default=50, help="The number of training data.")
    parser.add_argument('--M', type=int, default=3,  help="The number of hidden units")
    args = parser.parse_args()

    #=== データ数 ===
    N = args.N

    #=== ニューラルネットワークの形 ===
    D = 1       # 入力の次元(x座標) <簡単のため固定>
    M = args.M  # 隠れ層の数(基底関数の数)
    K = 1       # 出力層の数(y座表) <簡単のため固定>

    #=== 重みパラメータ ===
    W1  = M*(D+1) # 第1層の重みパラメータ数(ダミー変数を足す) # 隠れ層の重みは M*(D+1) 行列w1で表現
    W2  = K*(M+1) # 第2層の重みパラメータ数(ダミー変数を足す) # 出力層の重みは K*(M+1) 行列w2で表現
    W   = W1 + W2 # 重みパラメータ数

    #=== 最急降下法のパラメータ ===
    ALPHA = 0.01      # 最急降下法の勾配係数
    ITER_MAX = 5000   # 最大反復回数
    ITER_EPS = 1.0e-4 # 勾配のノルムがこれかになったら停止

    fit(args.func, args.lambda_func, args.xmin, args.xmax)

