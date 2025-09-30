\# Tutorial de Uso



A continuación se presentan ejemplos concretos que ilustran cómo emplear el módulo `gauss\_quadrature.py` para aproximar integrales definidas de funciones complejas, combinando la precisión de la cuadratura Gauss-Legendre con la elegancia del análisis simbólico.



```python

from gauss\_quadrature import funcion, obtener\_puntos\_pesos\_reescalados, funcion\_taylor

import numpy as np



\# Aproximar la integral de f(x) = x^6 - x^2 sin(2x) en el intervalo \[1, 3]

\# Usando 4 nodos de Gauss-Legendre para optimizar precisión con un mínimo de evaluaciones

puntos, pesos = obtener\_puntos\_pesos\_reescalados(1, 3, 4)

resultado = sum(pesos \* funcion(puntos))

print("Integral aproximada:", resultado)



\# Ahora, aproximamos la misma función mediante expansión de Taylor de grado 3

\# Esta técnica transforma el término trascendental en un polinomio, compatible con la cuadratura

valores, expr = funcion\_taylor(np.array(\[1.0, 2.0, 3.0]), 3)

print("Valores aproximados con Taylor:", valores)

print("Expresión simbólica de la función truncada:", expr)



