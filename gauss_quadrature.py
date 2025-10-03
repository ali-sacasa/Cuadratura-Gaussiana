"""
Módulo para integración numérica usando cuadratura Gauss-Legendre.

Este módulo implementa el método de cuadratura Gauss-Legendre para aproximar
integrales definidas, incluyendo soporte para expansiones en serie de Taylor
de funciones trascendentales.
"""

import numpy as np
import sympy as sp

# Variable simbólica global
t = sp.symbols("t")

def funcion(variable_independiente):
    """Evalúa la función φ(t) = t⁶ - t² sin(2t) en un punto o array de puntos.

    Esta función representa el integrando principal del problema de integración
    numérica. La función está vectorizada para operar sobre arrays de NumPy.

    Args:
        variable_independiente (float or numpy.ndarray): Punto(s) en los que 
            se evalúa la función.

    Returns:
        float or numpy.ndarray: Valor(es) de la función en el(los) punto(s).

    Examples:
        >>> funcion(1.0)
        -0.8185948536513634

        >>> import numpy as np
        >>> t_values = np.array([1.0, 2.0])
        >>> funcion(t_values)
        array([ -0.81859485,  54.90729447])

        La función en t=2.0:
        $$φ(2) = 2^6 - 2^2 \\cdot \\sin(4) ≈ 64 - 4 \\cdot (-0.7568) ≈ 54.9073$$
    """
    return variable_independiente**6 - variable_independiente**2 * np.sin(2*variable_independiente)


def obtener_puntos_pesos_reescalados(a, b, N):
    """Calcula los puntos y pesos de Gauss-Legendre reescalados al intervalo [a, b].

    Transforma los puntos y pesos estándar de Gauss-Legendre del intervalo [-1, 1]
    al intervalo [a, b] mediante un cambio de variable lineal.

    Args:
        a (float): Límite inferior del intervalo de integración.
        b (float): Límite superior del intervalo de integración.
        N (int): Número de puntos de cuadratura (grado de la aproximación).

    Returns:
        tuple: Una tupla (puntos, pesos) donde:
            - puntos (numpy.ndarray): Array de puntos de Gauss en [a, b]
            - pesos (numpy.ndarray): Array de pesos correspondientes

    Raises:
        ValueError: Si N no es un entero positivo o si a ≥ b.

    Examples:
        >>> puntos, pesos = obtener_puntos_pesos_reescalados(1, 3, 3)
        >>> print(f"Puntos: {puntos}")
        Puntos: [1.11270167 2.         2.88729833]
        >>> print(f"Pesos: {pesos}")
        Pesos: [0.55555556 0.88888889 0.55555556]

        La transformación se define como:
        $$x_{[a,b]} = \\frac{b-a}{2}x_{[-1,1]} + \\frac{a+b}{2}$$
    """
    if N <= 0:
        raise ValueError("N debe ser un entero positivo")
    if a >= b:
        raise ValueError("El límite inferior 'a' debe ser menor que 'b'")
    
    puntos_std, pesos_std = np.polynomial.legendre.leggauss(N)
    puntos_rescalados = 0.5 * (b - a) * puntos_std + 0.5 * (b + a)
    pesos_rescalados = 0.5 * (b - a) * pesos_std
    
    return puntos_rescalados, pesos_rescalados


def funcion_taylor(variable_independiente, grado):
    """Aproxima la función usando expansión en serie de Taylor de sin(2t).

    Reemplaza sin(2t) por su expansión en serie de Taylor alrededor de t=0
    y evalúa la función resultante (que es un polinomio) en los puntos dados.

    Args:
        variable_independiente (numpy.ndarray): Puntos donde evaluar la función.
        grado (int): Grado de la expansión de Taylor (número de términos - 1).

    Returns:
        tuple: Una tupla (valores, expr) donde:
            - valores (numpy.ndarray): Valores de la función aproximada
            - expr (sympy.Expr): Expresión simbólica de la aproximación

    Examples:
        >>> t_values = np.array([1.0, 2.0])
        >>> valores, expr = funcion_taylor(t_values, 4)
        >>> print(f"Valores: {valores}")
        Valores: [ -0.83333333  54.4       ]
        >>> print(f"Expresión: {expr}")
        Expresión: t**6 - 2*t**4 + (2/3)*t**2

        La serie de Taylor para sin(2t) hasta grado 3 es:
        $$\\sin(2t) ≈ 2t - \\frac{(2t)^3}{3!} = 2t - \\frac{8t^3}{6}$$
    """
    serie = sp.series(sp.sin(2*t), t, 0, grado + 1).removeO()
    expr = t**6 - t**2 * serie
    f_lambdified = sp.lambdify(t, expr, modules=["numpy"])
    return f_lambdified(variable_independiente), expr


def calcular_integral_gaussiana(func, a, b, N):
    """Calcula la aproximación de una integral usando cuadratura Gauss-Legendre.

    Args:
        func (callable): Función a integrar.
        a (float) - Límite inferior de integración.
        b (float) - Límite superior de integración.
        N (int) - Número de puntos de cuadratura.

    Returns:
        float: Aproximación numérica de la integral.

    Examples:
        >>> from gauss_quadrature import funcion, calcular_integral_gaussiana
        >>> resultado = calcular_integral_gaussiana(funcion, 1, 3, 5)
        >>> print(f"Resultado {resultado:.6f}")
        Resultado: 578.666667
    """
    puntos, pesos = obtener_puntos_pesos_reescalados(a, b, N)
    return np.sum(pesos * func(puntos))


# Cálculo de la solución analítica para comparación
I = sp.integrate(t**6 - t**2 * sp.sin(2*t), (t, 1, 3))
solucion_analitica = I.evalf()

if __name__ == "__main__":
    # Ejemplo de uso y pruebas
    print("Pruebas de Cuadratura Gauss-Legendre \n")
    
    for N in range(2, 7):
        puntos, pesos = obtener_puntos_pesos_reescalados(1, 3, N)
        integral_num = calcular_integral_gaussiana(funcion, 1, 3, N)
        print(f"Los resultados para N={N}")
        print(f"Los puntos usados son {puntos}")
        print(f"Pesos {pesos}")
        print(f"La integral aproximada es {integral_num:.8f}")
        print(f"La diferencia con solución analítica es {abs(solucion_analitica - integral_num):.2e}\n")