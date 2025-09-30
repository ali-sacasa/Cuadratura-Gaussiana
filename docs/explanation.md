\# Método de Cuadratura Gauss-Legendre



En física y matemáticas aplicadas, a menudo necesitamos evaluar integrales definidas de funciones continuas cuya forma explícita impide una solución analítica directa.  

La \*\*cuadratura Gauss-Legendre\*\* ofrece un método elegante y eficiente para aproximar estas integrales:



$$

\\int\_a^b f(x) \\, dx \\approx \\sum\_{i=1}^{N} w\_i \\, f(x\_i),

$$



donde los $x\_i$ denominados \*\*nodos de Gauss\*\*, se seleccionan cuidadosamente dentro del intervalo de integración $\[a, b]$, y los $w\_i$, llamados \*\*pesos de Gauss\*\*, se asignan de manera que la aproximación sea exacta para polinomios de grado hasta $2N-1$.  



Este enfoque no solo minimiza el número de evaluaciones de la función, sino que también garantiza una alta precisión en aplicaciones físicas, por ejemplo, al calcular momentos de inercia, energía potencial en campos continuos o integrales de propagadores en mecánica cuántica.  



Para funciones que involucran términos trascendentales, como $\\sin(2x)$, podemos aplicar la \*\*expansión en serie de Taylor\*\*, convirtiendo la función en un polinomio truncado. Posteriormente, la integral polinómica se evalúa usando los nodos y pesos de Gauss-Legendre, combinando el análisis simbólico con la eficiencia numérica de la cuadratura.



