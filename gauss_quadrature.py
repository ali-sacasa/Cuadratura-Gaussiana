# Sebastián Alí Sacasa Céspedes - C4J546
# Tarea II - Cuadratura de Gauss-Legendre
# Fecha: 28 de octubre del 2023
# Correo del autor correspondiente: sebastian.sacasa@ucr.ac.cr

# Se importan las librerías necesarias para los cálculos numéricos y simbólicos.
import numpy as np
import sympy as sp
# Se define la función a integrar con el uso de operaciones de numpy asegura la vectorización. 
# Las variables se redactan en formato snake_case.

# Cálculo simbólico de la integral para comparar resultados
# Se usa sympy para definir la variable simbólica y calcular la integral definida.
t = sp.symbols("t")
I = sp.integrate(t**6 - t**2*sp.sin(2*t), (t, 1, 3))
solucion_analítica = I.evalf()

def funcion(variable_independiente):
    return variable_independiente**6 -variable_independiente**2*np.sin(2*variable_independiente)

#Esta función toma el número de puntos N y los límites de integración [a, b] y devuelve los puntos y pesos apropiados para ese intervalo.
# Básicamente es una transformación lineal de los puntos y pesos estándar de Gauss-Legendre en [-1, 1] al intervalo [a, b], un cambio de variable.
def obtener_puntos_pesos_reescalados(a, b, N):

    # np.polynomial.legendre.leggauss(N) devuelve los puntos y pesos en el intervalo estándar [-1, 1].
    puntos_std, pesos_std = np.polynomial.legendre.leggauss(N)
    
    # Reescalado de los puntos y pesos al intervalo [a, b]
    puntos_rescalados = 0.5*(b-a)*puntos_std+0.5*(b+a)
    
    pesos_rescalados =0.5*(b-a)*pesos_std
    
    return puntos_rescalados, pesos_rescalados


# Se definen los límites de integración, en este caso, no se realizará la integral en el intervalo [-1, 1], sino en [0, 2].
# A su vez, no se analizarán las desigualdades para los límites de integración, ya en la función de obtener_puntos_pesos_reescalados(a, b, N),
# al ser lineal, todo cambio se reduce a una transformación lineal o cambio de orden integración.


for N in range(2, 7):
    puntos, pesos = obtener_puntos_pesos_reescalados(1, 3, N)
    integral_num = np.sum(pesos * funcion(puntos))
    print(f"Resultados para N={N}")
    print(f"Puntos de colocación = {puntos}")
    print(f"Pesos = {pesos}")
    print(f"Resultado de la integral = {integral_num}\n")

print("Comparación con la solución analítica")
for N in range(2, 7):
    puntos, pesos=obtener_puntos_pesos_reescalados(1, 3, N)
    integral_num=np.sum(pesos*funcion(puntos))
    print(f"Diferencia con N={N}: {abs(solucion_analítica-integral_num)}")


# Expansión de sin(2t) en serie de Taylor y búsqueda de N óptimo con error menor o igual a 1e-18
error_tolerancia = 1e-18
def funcion_taylor(variable_independiente, grado):
    serie = sp.series(sp.sin(2*t), t, 0, grado+1).removeO()
    expr = t**6 - t**2 * serie
    f_lambdified = sp.lambdify(t, expr, modules=["numpy"])
    return f_lambdified(variable_independiente), expr
grado = 2
encontrado = False
maximo_grado = 50  # límite para evitar bucles infinitos en caso de no encontrar solución
while grado<= maximo_grado:  # bucle infinito, salimos cuando encontremos la condición
    valores_funcion, expr_taylor = funcion_taylor(np.array([0.0]), grado)  # solo para construir expr
    integral_analitica_taylor = sp.integrate(expr_taylor, (t, 1, 3)).evalf()

    for N in range(2, grado+3):  # N debe crecer con el grado
        puntos, pesos = obtener_puntos_pesos_reescalados(a=1, b=3, N=N)
        valores, _ = funcion_taylor(puntos, grado)
        integral_num = np.sum(pesos * valores)
        error = float(abs(integral_analitica_taylor - integral_num))
        if error <= error_tolerancia:
            print(f"Grado = {grado}, N = {N}, error = {error:.2e} <= {error_tolerancia}")
            print(f"Integral analítica (polinómica) = {integral_analitica_taylor}")
            print(f"Integral numérica (Gauss-Leg)  = {integral_num}\n")
            encontrado = True
            break
    if encontrado:
        break
    grado += 1
if not encontrado:
    print(f"No se encontró N para grado hasta {maximo_grado} con error <= {error_tolerancia}")
#El valor óptimo de N depende del grado de la expansión en serie de Taylor utilizada para aproximar sin(2t).
#A medida que se incrementa el grado de la serie de Taylor, se requiere un mayor número de puntos N en la cuadratura de Gauss-Legendre para alcanzar el mismo nivel de precisión.
#Se ajusta una cota para evitar bucles infinitos en caso de no encontrar una solución dentro de un rango razonable que ocupe mucha memoria y tiempo de procesamiento.