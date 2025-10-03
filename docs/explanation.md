# Método de Cuadratura Gauss-Legendre

En física y matemáticas aplicadas, a menudo necesitamos evaluar **integrales definidas** de funciones continuas cuya forma explícita impide una solución analítica directa. Este tipo de problemas aparece, por ejemplo, al calcular momentos de inercia, energías potenciales en campos continuos o integrales de propagadores en mecánica cuántica.  

La **cuadratura Gauss-Legendre** ofrece un método elegante y eficiente para aproximar estas integrales

$$
\int_a^b \phi(x) \, dx \approx \sum_{i=1}^{N} w_i \, \phi(x_i)
$$

donde los \(x_i\), denominados **nodos de Gauss**, se seleccionan cuidadosamente dentro del intervalo \([a,b]\), y los \(w_i\), llamados **pesos de Gauss**, se asignan de manera que la aproximación sea exacta para **polinomios de grado hasta \(2N-1\)**. Esto permite obtener resultados muy precisos incluso con un número reducido de nodos.  

---

## Aplicación a funciones trascendentales

Cuando la función a integrar involucra términos trascendentales, como \(\sin(2x)\), podemos recurrir a la **expansión en serie de Taylor**, truncando la serie a un polinomio de grado adecuado. La integral resultante, ahora polinómica, se evalúa usando los nodos y pesos de Gauss-Legendre, combinando el rigor del análisis simbólico con la eficiencia numérica de la cuadratura.  

Sea por ejemplo

$$
\phi(t) = t^6 - t^2 \sin(2t)
$$

y consideremos su aproximación mediante una serie de Taylor de grado \(k\)

$$
\sin(2t) \approx \sum_{n=0}^{k} \frac{(-1)^n (2t)^{2n+1}}{(2n+1)!}.
$$

La integral se convierte entonces en un polinomio que puede evaluarse **exactamente** con la cuadratura de Gauss siempre que el grado del polinomio cumpla \(G[\phi(t)] \leq 2N-1\).

---

## Selección del número de nodos \(N\)

La elección de \(N\) es crítica para asegurar la exactitud numérica

- Para \(N=3\), la regla de Gauss-Legendre es exacta para polinomios de grado hasta \(2\cdot 3-1 = 5\).  
- Dado que el polinomio resultante de la expansión de Taylor tiene grado 4, se cumple \(4 \leq 5\), y la integral se calcula exactamente dentro de la **precisión de máquina**.  
- Para \(N=2\), el grado máximo exacto es \(2\cdot 2-1 = 3\), por lo que \(4 \not\leq 3\), y la integral es solo aproximada, con un error de orden \(\mathcal{O}(10^{-1})\).  

Esto muestra cómo la **regla de precisión de la cuadratura** conecta la teoría de polinomios de Legendre con la eficiencia numérica.

---

## Ejemplo práctico en Python

A continuación, se muestra un ejemplo de cómo usar el módulo `gauss_quadrature.py` para aproximar la integral de \(\phi(t) = t^6 - t^2 \sin(2t)\) en el intervalo \([1,3]\):

"Python"

import numpy as np
from gauss_quadrature import funcion, obtener_puntos_pesos_reescalados, funcion_taylor

"Número de nodos"
N = 4

Puntos y pesos reescalados
puntos, pesos = obtener_puntos_pesos_reescalados(1, 3, N)

"Integral usando la función original"

resultado = np.sum(pesos * funcion(puntos))
print("Integral aproximada", resultado)

"Integral usando expansión de Taylor de grado 3"

valores, expr = funcion_taylor(puntos, 3)
resultado_taylor = np.sum(pesos * valores)
print("Integral aproximada con Taylor", resultado_taylor)
print("Expresión simbólica truncada", expr)
