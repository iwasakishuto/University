#coding: utf-8
import argparse
import numpy as np
import matplotlib.pyplot as plt
from sklearn.linear_model import LinearRegression
from sklearn.metrics import mean_squared_error as mse

if __name__ == "__main__":
  tp = lambda x:list(map(float, x.split(',')))

  parser = argparse.ArgumentParser()
  parser.add_argument('--Concentration', type=tp,  help='Input the concentration, with ",". ex.) --Concentration  0,0.10,...,0.50')
  parser.add_argument('--Absorbance',    type=tp,  help='Input the result absorbance, with ",". ex.) --Absorbance  0.048,0.045,...,0.076')
  parser.add_argument('--unit',          type=str, help='Specify the unit of the Concentration. ex.) mg/mL', default='mg/mL')
  parser.add_argument('--wavelength',    type=str, help='Specify the wavelength of the Absorbance. ex) 540')
  args = parser.parse_args()

  c = np.array(args.Concentration).reshape(-1, 1)
  a = np.array(args.Absorbance).reshape(-1, 1)
  unit = args.unit
  wavelength = args.wavelength

  lr = LinearRegression()
  lr.fit(c,a)

  coef = lr.coef_[0][0]
  intercept = lr.intercept_[0]

  Calibration_curve = "$y={:.3f}x{}{:.3f}$".format(coef, "+" if intercept > 0 else "",  intercept) # 検量線
  Coefficient_of_determination = "$R^2$ = {:.3f}".format(lr.score(c,a)) # 決定係数
  Root_Mean_Square_Error = "RMSE = {:.3f}".format(mse(a, lr.predict(c)) ** (1/2)) # 平均二乗誤差

  plt.plot(c, lr.predict(c), label="Calibration curve: "+Calibration_curve)
  plt.scatter(c,a, label="{}\n{}".format(Coefficient_of_determination,
                                        Root_Mean_Square_Error))
  plt.xlabel("Concentration [{}]".format(unit))
  plt.ylabel("Absorbance ({})".format(wavelength))
  plt.title("The relationship between Concentration and Absorbance")
  plt.legend()
  plt.savefig("calibration-curve.png")

  lr_inv = LinearRegression()
  lr_inv.fit(a,c)
  coef_inv = lr_inv.coef_[0][0]
  intercept_inv = lr_inv.intercept_[0]
  a2c = "y={:.3f}x{}{:.3f}".format(coef_inv, "+" if intercept_inv > 0 else "",  intercept_inv)
  print("Clibration curve (Absorbance → Concentration): {}".format(a2c))
