# Tutorial de Uso – Cuadratura Gauss-Legendre

En este tutorial se muestra cómo utilizar el módulo "gauss_quadrature.py" para aproximar integrales mediante la **cuadratura Gauss-Legendre**, incluyendo la aplicación de **expansión de Taylor** para funciones trascendentales.

---

## Aproximación de una integral definida

Consideremos la función:

$$
f(t) = t^6 - t^2 \sin(2t),
$$

y el intervalo de integración \([1,3]\). Podemos aproximar su integral usando \(N\) nodos de Gauss-Legendre:

```python
from gauss_quadrature import funcion, obtener_puntos_pesos_reescalados
import numpy as np

# Número de nodos
N = 4

# Obtener puntos y pesos reescalados al intervalo [1,3]
puntos, pesos = obtener_puntos_pesos_reescalados(1, 3, N)

# Evaluar la integral numéricamente
resultado = np.sum(pesos * funcion(puntos))
print("Integral aproximada:", resultado)
