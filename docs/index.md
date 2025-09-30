# Introducción

En física y matemáticas aplicadas, con frecuencia nos enfrentamos a la necesidad de evaluar **integrales definidas** de funciones continuas cuyo cálculo analítico resulta complejo o incluso imposible de realizar de forma directa. Estas integrales aparecen en múltiples contextos:  

- Cálculo de momentos de inercia en cuerpos continuos.  
- Evaluación de energías potenciales en campos físicos.  
- Integración de propagadores y funciones de Green en mecánica cuántica.  

Para abordar este tipo de problemas, es fundamental combinar **rigor matemático** con **eficiencia computacional**. Una de las herramientas más poderosas para este propósito es la **cuadratura de Gauss-Legendre**, que permite aproximar integrales definidas con una **precisión controlada** usando un número limitado de evaluaciones de la función.  

En esta documentación se explorará:  

1. El **método de cuadratura Gauss-Legendre** y su fundamentación teórica.  
2. La **aplicación a funciones trascendentales** mediante expansiones en serie de Taylor.  
3. Ejemplos prácticos de uso en Python utilizando el módulo `gauss_quadrature.py`.  
4. Una **referencia completa de las funciones** implementadas y sus docstrings.  

El objetivo de este proyecto es proporcionar una guía clara y completa, que combine teoría, implementación y ejemplos, para que el lector pueda comprender y aplicar la **cuadratura Gauss-Legendre** en contextos de física aplicada y matemáticas computacionales.Esto para efectos del curso de Física Computacional FS0432, específicamente la tarea II y otros fines más allá de este, puesto que, la integración numérica es una herramienta fundamental en física y matemáticas aplicadas.  

