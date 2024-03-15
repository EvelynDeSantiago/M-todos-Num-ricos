# Métodos Numéricos
Los métodos numéricos son técnicas utilizadas para resolver problemas matemáticos mediante aproximaciones numéricas en lugar de soluciones analíticas exactas. Estos métodos son especialmente útiles cuando las soluciones analíticas no son factibles o son muy complejas.
## Introducción
Los métodos numéricos constituyen una herramienta esencial en la resolución de problemas matemáticos y científicos en los cuales las soluciones analíticas exactas son difíciles o imposibles de obtener. Estos métodos se basan en la aproximación de soluciones mediante técnicas computacionales, permitiendo obtener resultados aceptables con un grado controlado de error.

En este contexto, el objetivo fundamental de los métodos numéricos es desarrollar algoritmos eficientes y precisos para resolver problemas matemáticos, tales como la integración numérica, la solución de ecuaciones diferenciales, la interpolación de datos, entre otros. Estos problemas son comunes en diversas áreas de la ciencia, la ingeniería, la economía y muchas otras disciplinas.

En esta introducción, exploraremos algunos conceptos fundamentales de los métodos numéricos, incluyendo la representación de números en la computadora, el análisis de errores, y las principales técnicas algorítmicas utilizadas en la resolución de problemas numéricos. Además, discutiremos la importancia de la convergencia y la estabilidad de los métodos numéricos, así como su aplicación en la modelización y simulación de fenómenos reales.

### TEMA 3
# Método Gauss-seidel
El método de Gauss-Seidel es un procedimiento iterativo para hallar soluciones aproximadas a un sistema de ecuaciones algebraicas lineales con una precisión arbitrariamente elegida. El método se aplica a matrices cuadradas con elementos no nulos en sus diagonales y la convergencia se garantiza si la matriz es diagonalmente dominante.
## Algoritmo:
Dado el sistema de ecuaciones:
9X1+2X2−X3=−2
7X1+8X2+5X3=3
3X1+4X2−10X3=6
1-.Definimos la matriz de coeficientes A y el vector de términos independientes B.
2-. Vamos a utilizar una aproximación inicial de X(0)=(000).
3-. Luego, aplicamos el método de Gauss-Seidel hasta que las soluciones converjan.
## Algoritmo Gauss_Seidel

		Definir A[3, 3] como Entero
		Definir B[3] como Entero
		Definir X[3] como Real
		Definir X_ant[3] como Real
		Definir tol como Real
		Definir iterMax como Entero
		Definir n como Entero
		Definir iter como Entero
		Definir i como Entero
		Definir j como Entero
		Definir suma como Real
		
		// Inicializar la matriz A y el vector B
		A[1][1] <- 9
		A[1][2] <- 2
		A[1][3] <- -1
		A[2][1] <- 7
		A[2][2] <- 8
		A[2][3] <- 5
		A[3][1] <- 3
		A[3][2] <- 4
		A[3][3] <- -10
		B[1] <- -2
		B[2] <- 3
		B[3] <- 6
		
		// Especificar la tolerancia y el número máximo de iteraciones
		tol <- 0.0001
		iterMax <- 100
		n <- Longitud(A)
		
		// Inicializar soluciones
		Para i <- 1 Hasta n
			X[i] <- 0
		FinPara
		
		iter <- 0
		
		// Bucle principal del método de Gauss-Seidel
		Mientras (iter < iterMax) Y (No Converge(X, X_ant, tol))
			// Copiar las soluciones de la iteración anterior
			Para i <- 1 Hasta n
				X_ant[i] <- X[i]
			FinPara
			
			// Actualizar cada solución
			Para i <- 1 Hasta n
				suma <- 0
				Para j <- 1 Hasta n
					Si (j <> i) Entonces
						suma <- suma + A[i][j] * X[j]
					FinSi
				FinPara
				X[i] <- (B[i] - suma) / A[i][i]
			FinPara
			
			iter <- iter + 1
		FinMientras
		
		// Mostrar las soluciones
		Escribir "Soluciones:"
		Para i <- 1 Hasta n
			Escribir "X" + i + ": " + X[i]
		FinPara
FinSubproceso

Funcion Converge(X[] como Real, X_ant[] como Real, tol como Real) como Logico
    Definir n como Entero
    Definir i como Entero
    n <- Longitud(X)
    Para i <- 1 Hasta n
        Si (Absoluto(X[i] - X_ant[i]) > tol) Entonces
            Devolver Falso
        FinSi
    FinPara
    Devolver Verdadero
FinFuncion

## Implementación de usos:
Una vez que tienes los datos organizados, debes crear las fórmulas necesarias para implementar el método Gauss-Seidel. En cada celda correspondiente a una variable, debes ingresar una fórmula que tome en cuenta las demás variables y sus respectivos valores. Para ello, puedes utilizar las funciones de Excel, como SUMA, RESTA, PRODUCTO, entre otras.

## Ejercicio 1:
package Ejercicios;

import java.util.Arrays;

public class GaussSeidel {

    public static double[] gaussSeidel(double[][] coeficientes, double[] terminosIndependientes, double tol, int iterMax) {
        int n = terminosIndependientes.length;
        double[] soluciones = new double[n];
        double[] solucionesAnterior = new double[n];
        int iter = 0;

        // Inicio de las soluciones
        for (int i = 0; i < n; i++) {
            soluciones[i] = 0.0;
        }

        // Bucle principal del método de Gauss-Seidel
        do {
            // Copiar las soluciones de la iteración anterior
            System.arraycopy(soluciones, 0, solucionesAnterior, 0, n);

            // Actualizar cada solución
            for (int i = 0; i < n; i++) {
                double suma = 0.0;
                for (int j = 0; j < n; j++) {
                    if (j != i) {
                        suma += coeficientes[i][j] * soluciones[j];
                    }
                }
                soluciones[i] = (terminosIndependientes[i] - suma) / coeficientes[i][i];
            }

            iter++;
        } while (iter < iterMax && !converge(soluciones, solucionesAnterior, tol));

        return soluciones;
    }

    // Función para verificar si las soluciones han convergido
    public static boolean converge(double[] soluciones, double[] solucionesAnterior, double tol) {
        for (int i = 0; i < soluciones.length; i++) {
            if (Math.abs(soluciones[i] - solucionesAnterior[i]) > tol) {
                return false;
            }
        }
        return true;
    }

    public static void main(String[] args) {
        // Definir la matriz de coeficientes
        double[][] coeficientes = {
                {9, 2, -1},
                {7, 8, 5},
                {3, 4, -10}
        };

        // Definir el vector de términos independientes
        double[] terminosIndependientes = {-2, 3, 6};

        // Especificar la tolerancia y el número máximo de iteraciones
        double tol = 0.0001;
        int iterMax = 100;

        // Resolver el sistema de ecuaciones lineales
        double[] soluciones = gaussSeidel(coeficientes, terminosIndependientes, tol, iterMax);

        // Impreción de las soluciones
        System.out.println("Soluciones:");
        for (int i = 0; i < soluciones.length; i++) {
            System.out.println("x" + (i + 1) + ": " + soluciones[i]);
        }
    }
}
## Ejercicio 2:
public class GaussSeidel {
    // Función para imprimir la matriz
    public static void imprimirMatriz(double[][] matriz) {
        int n = matriz.length;
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n + 1; j++) {
                System.out.print(matriz[i][j] + "\t");
            }
            System.out.println();
        }
    }

    public static double[] gaussSeidel(double[][] matriz, double tolerancia, int iterMax) {
        int n = matriz.length;
        double[] x = new double[n];
        double[] xAnterior = new double[n];
        int iter = 0;

 
        for (int i = 0; i < n; i++) {
            x[i] = 0;
        }

       
        while (iter < iterMax) {
            
            for (int i = 0; i < n; i++) {
                xAnterior[i] = x[i];
            }

            
            for (int i = 0; i < n; i++) {
                double sum = 0;
                for (int j = 0; j < n; j++) {
                    if (j != i) {
                        sum += matriz[i][j] * x[j];
                    }
                }
                x[i] = (matriz[i][n] - sum) / matriz[i][i];
            }

            
            double error = 0;
            for (int i = 0; i < n; i++) {
                error += Math.abs(x[i] - xAnterior[i]);
            }

            if (error < tolerancia) {
                System.out.println("Convergencia alcanzada en la iteración " + (iter + 1));
                return x;
            }

            iter++;
        }

        System.out.println("Iteraciones máximas alcanzadas sin convergencia.");
        return null;
    }

    public static void main(String[] args) {
        // Ejemplo de sistema de ecuaciones:
        // 3x + y - z = 5
        // x - 4y + z = -7
        // 2x + y - 3z = 3
        double[][] matriz = {
            {3, 1, -1, 5},
            {1, -4, 1, -7},
            {2, 1, -3, 3}
        };

        // Tolerancia y número máximo de iteraciones
        double tolerancia = 0.0001;
        int iterMax = 100;

        // Llamar al método de Gauss-Seidel
        double[] soluciones = gaussSeidel(matriz, tolerancia, iterMax);

        // Imprimir las soluciones
        if (soluciones != null) {
            System.out.println("Soluciones:");
            for (int i = 0; i < soluciones.length; i++) {
                System.out.println("x[" + i + "] = " + soluciones[i]);
            }
        }
    }
}

## Ejercicio 3:
public class GaussSeidel {
    
    static final int MAX_ITERATIONS = 100;
    static final double TOLERANCE = 0.0001;

    public static void main(String[] args) {
        double[][] coefficients = {
            {10, -1, 2, 0},
            {-1, 11, -1, 3},
            {2, -1, 10, -1},
            {0, 3, -1, 8}
        };
        double[] initialValues = {0, 0, 0, 0};
        double[] solutions = gaussSeidel(coefficients, initialValues);
        System.out.println("Solutions:");
        for (int i = 0; i < solutions.length; i++) {
            System.out.println("x[" + i + "] = " + solutions[i]);
        }
    }

    public static double[] gaussSeidel(double[][] coefficients, double[] initialValues) {
        int n = coefficients.length;
        double[] solutions = new double[n];
        double[] newSolutions = new double[n];
        int iterations = 0;
        double error;

        do {
            for (int i = 0; i < n; i++) {
                double sum = 0.0;
                for (int j = 0; j < n; j++) {
                    if (j != i) {
                        sum += coefficients[i][j] * newSolutions[j];
                    }
                }
                newSolutions[i] = (initialValues[i] - sum) / coefficients[i][i];
            }

            error = 0.0;
            for (int i = 0; i < n; i++) {
                error = Math.max(error, Math.abs(newSolutions[i] - solutions[i]));
                solutions[i] = newSolutions[i];
            }

            iterations++;
        } while (error > TOLERANCE && iterations < MAX_ITERATIONS);

        if (iterations >= MAX_ITERATIONS) {
            System.out.println("Maximum iterations reached without convergence");
        }

        return solutions;
    }
}
## Ejercicio 4:
public class GaussSeidel {
    
    public static void main(String[] args) {
        double[][] coeficientes = {{12, 3, -5}, {3, -8, -3}, {5, 4, -12}}; // Coeficientes de las variables
        double[] constantes = {1, 1, 3}; // Términos constantes del sistema de ecuaciones
        double[] solucionInicial = {0, 0, 0}; // Suposición inicial para las soluciones
        
        int iteracionesMaximas = 100; // Número máximo de iteraciones
        double tolerancia = 0.0001; // Tolerancia para la convergencia
        
        double[] solucion = gaussSeidel(coeficientes, constantes, solucionInicial, iteracionesMaximas, tolerancia);
        
        // Imprimir la solución
        System.out.println("Solución:");
        for (int i = 0; i < solucion.length; i++) {
            System.out.println("x" + (i+1) + " = " + solucion[i]);
        }
    }
    
    // Método de Gauss-Seidel para resolver sistemas de ecuaciones lineales
    public static double[] gaussSeidel(double[][] coeficientes, double[] constantes, double[] solucionInicial, int iteracionesMaximas, double tolerancia) {
        int n = constantes.length;
        double[] solucion = new double[n];
        double[] solucionAnterior = new double[n];
        int iteraciones = 0;
        double error = Double.MAX_VALUE;
        
        // Copiar la solución inicial
        for (int i = 0; i < n; i++) {
            solucion[i] = solucionInicial[i];
        }
        
        // Iterar hasta que se alcance la tolerancia o el número máximo de iteraciones
        while (error > tolerancia && iteraciones < iteracionesMaximas) {
            for (int i = 0; i < n; i++) {
                solucionAnterior[i] = solucion[i];
            }
            
            for (int i = 0; i < n; i++) {
                double suma = 0.0;
                for (int j = 0; j < n; j++) {
                    if (j != i) {
                        suma += coeficientes[i][j] * solucion[j];
                    }
                }
                solucion[i] = (constantes[i] - suma) / coeficientes[i][i];
            }
            
            // Calcular el error
            error = 0.0;
            for (int i = 0; i < n; i++) {
                error = Math.max(error, Math.abs(solucion[i] - solucionAnterior[i]));
            }
            
            iteraciones++;
        }
        
        if (error <= tolerancia) {
            System.out.println("Convergencia alcanzada en " + iteraciones + " iteraciones.");
        } else {
            System.out.println("El método no converge en " + iteracionesMaximas + " iteraciones.");
        }
        
        return solucion;
    }
}
## Ejercicio 5:
public class GaussSeidelEjemplo5 {

    public static void main(String[] args) {
        double[][] sistemaEcuaciones = {
                {5, 1, -2, 19},
                {-3, -4, 1, -28},
                {2, -1, 6, 31}
        };

        double[] solucionesIniciales = {0, 0, 0};

        resolverSistema(sistemaEcuaciones, solucionesIniciales, 0.0001, 50);
    }

    public static void resolverSistema(double[][] matriz, double[] soluciones, double tolerancia, int iteracionesMax) {
        int filas = matriz.length;

        for (int iteracion = 1; iteracion <= iteracionesMax; iteracion++) {
            double[] nuevasSoluciones = new double[filas];

            for (int i = 0; i < filas; i++) {
                double suma = 0;
                for (int j = 0; j < filas; j++) {
                    if (j != i) {
                        suma += matriz[i][j] * soluciones[j];
                    }
                }
                nuevasSoluciones[i] = (matriz[i][filas] - suma) / matriz[i][i];
            }

            // Comprobar convergencia
            if (convergencia(soluciones, nuevasSoluciones, tolerancia)) {
                System.out.println("Soluciones convergen en la iteración " + iteracion + ":");
                mostrarSoluciones(nuevasSoluciones);
                return;
            }

            soluciones = nuevasSoluciones;
        }

        System.out.println("No se alcanzó convergencia después de " + iteracionesMax + " iteraciones.");
    }

    private static boolean convergencia(double[] anteriores, double[] actuales, double tolerancia) {
        for (int i = 0; i < anteriores.length; i++) {
            if (Math.abs(actuales[i] - anteriores[i]) > tolerancia) {
                return false;
            }
        }
        return true;
    }

    private static void mostrarSoluciones(double[] soluciones) {
        for (int i = 0; i < soluciones.length; i++) {
            System.out.println("x" + (i + 1) + " = " + soluciones[i]);
        }
    }
}

# Método de Jacobi
En análisis numérico el método de Jacobi es un método iterativo, usado para resolver sistemas de ecuaciones lineales del tipo . El algoritmo toma su nombre del matemático alemán Carl Gustav Jakob Jacobi. El método de Jacobi consiste en usar fórmulas como iteración de punto fijo.
## ALGORITMO
6x+2y− z=4   x+5y+ z=3    2x+ y+4z=27 
1-.Inicializar valores iniciales para X, Y y Z.
2-.Para cada ecuación i, calcular xi utilizando los valores actuales de xj, donde j=i.
3-.Repetir el paso 2 hasta que se alcance la convergencia o se cumpla un número máximo de iteraciones.
## Algoritmo Gauss_Jacobi;
		// Definir variables
		Definir coeficientes[3, 3] como Entero
		Definir terminosIndependientes[3] como Entero
		Definir soluciones[3] como Real
		Definir solucionesAnterior[3] como Real
		Definir tol como Real
		Definir iterMax como Entero
		Definir iter como Entero
		
		// Inicializar matriz de coeficientes
		coeficientes[1][1] = 6
		coeficientes[1][2] = 2
		coeficientes[1][3] = -1
		coeficientes[2][1] = 1
		coeficientes[2][2] = 5
		coeficientes[2][3] = 1
		coeficientes[3][1] = 2
		coeficientes[3][2] = 1
		coeficientes[3][3] = 4
		
		// Inicializar vector de términos independientes
		terminosIndependientes[1] = 4
		terminosIndependientes[2] = 3
		terminosIndependientes[3] = 27
		
		// Especificar la tolerancia y el número máximo de iteraciones
		tol = 0.0001
		iterMax = 100
		iter = 0
		
		// Inicialización de las soluciones
		Para i desde 1 hasta 3
			soluciones[i] = 0.0
		FinPara
		
		// Bucle principal del método de Jacobi
		Mientras (iter < iterMax)
			// Copiar las soluciones de la iteración anterior
			Para i desde 1 hasta 3
				solucionesAnterior[i] = soluciones[i]
			FinPara
			
			// Actualizar cada solución
			Para i desde 1 hasta 3
				Real suma = 0.0
				Para j desde 1 hasta 3
					Si (j <> i) Entonces
						suma = suma + coeficientes[i][j] * solucionesAnterior[j]
					FinSi
				FinPara
				soluciones[i] = (terminosIndependientes[i] - suma) / coeficientes[i][i]
			FinPara
			
			iter = iter + 1
		FinMientras
		
		// Imprimir las soluciones
		Escribir("Soluciones:")
		Para i desde 1 hasta 3
			Escribir("x" + i + ": " + soluciones[i])
		FinPara
FinProceso


## IMPLEMENTACIÓN DE USO
Ahora implementamos el método de Jacobi para resolver un sistema de ecuaciones 
lineales  de  n  x  n  aproximando  la  solución  mediante  varias  iteraciones  sucesivas.  El 
lenguaje de programación utilizado será Java.

## EJERCICIO 1:
package Ejercicios;

public class MetodoJacobi {

    // Función para calcular el siguiente valor de x[i]
    public static double calculateNewX(double x, double y, double z) {
        return (4 - 2 * y + z) / 6.0;
    }

    // Función para calcular el siguiente valor de y[i]
    public static double calculateNewY(double x, double y, double z) {
        return (3 - x - z) / 5.0;
    }

    // Función para calcular el siguiente valor de z[i]
    public static double calculateNewZ(double x, double y, double z) {
        return (27 - 2 * x - y) / 4.0;
    }

    // Método de Jacobi
    public static void jacobi(double x, double y, double z, double tolerance, int maxIterations) {
        int iterations = 0;

        while (iterations < maxIterations) {
            double newX = calculateNewX(x, y, z);
            double newY = calculateNewY(x, y, z);
            double newZ = calculateNewZ(x, y, z);

            // Comprobar la convergencia
            double maxDiff = Math.max(Math.abs(newX - x), Math.max(Math.abs(newY - y), Math.abs(newZ - z)));

            if (maxDiff < tolerance) {
                System.out.println("Convergencia alcanzada en la iteración " + iterations);
                System.out.println("x = " + newX);
                System.out.println("y = " + newY);
                System.out.println("z = " + newZ);
                return;
            }

            // Actualizar los valores de x, y, y z
            x = newX;
            y = newY;
            z = newZ;

            iterations++;
        }

        System.out.println("No se alcanzó la convergencia después de " + maxIterations + " iteraciones.");
    }

    public static void main(String[] args) {
        double x = 0; // Valor inicial de x
        double y = 0; // Valor inicial de y
        double z = 0; // Valor inicial de z
        double tolerance = 0.0001; // Tolerancia
        int maxIterations = 100; // Número máximo de iteraciones

        jacobi(x, y, z, tolerance, maxIterations);
    }
}
## Ejercicio 2:
public class JacobiMethod {

    public static double[] jacobi(double[][] A, double[] b, int maxIterations, double tolerance) {
        int n = A.length;
        double[] x = new double[n];
        double[] x_new = new double[n];
        int iter = 0;
        double error = tolerance + 1; // Inicializar el error con un valor mayor que la tolerancia
        
        // Comenzar iteraciones
        while (iter < maxIterations && error > tolerance) {
            for (int i = 0; i < n; i++) {
                double sum = 0;
                for (int j = 0; j < n; j++) {
                    if (j != i) {
                        sum += A[i][j] * x[j];
                    }
                }
                x_new[i] = (b[i] - sum) / A[i][i];
            }
            
            // Calcular el error
            error = calculateError(x, x_new);
            
            // Actualizar el vector x con los nuevos valores
            System.arraycopy(x_new, 0, x, 0, n);
            
            iter++;
        }
        
        return x;
    }
    
    // Método para calcular el error entre dos iteraciones consecutivas
    private static double calculateError(double[] x, double[] x_new) {
        double max = Math.abs(x_new[0] - x[0]);
        for (int i = 1; i < x.length; i++) {
            double error = Math.abs(x_new[i] - x[i]);
            if (error > max) {
                max = error;
            }
        }
        return max;
    }

    public static void main(String[] args) {
        double[][] A = {{4, 1, 2}, {3, 5, 1}, {1, 1, 3}};
        double[] b = {4, 7, 3};
        int maxIterations = 1000;
        double tolerance = 1e-6;

        double[] solution = jacobi(A, b, maxIterations, tolerance);

        // Mostrar la solución
        System.out.println("Solución encontrada:");
        for (double value : solution) {
            System.out.println(value);
        }
    }
}
## Ejercicio 3:
import java.util.Arrays;

public class Jacobi {
    
    public static void main(String[] args) {
        double[][] coeficientes = {
            {10, -1, 2, 0},
            {-1, 11, -1, 3},
            {2, -1, 10, -1}
        };
        
        double[] soluciones = {6, 25, -11};
        
        double[] resultado = jacobi(coeficientes, soluciones, 0.0001, 100);
        
        System.out.println("Soluciones:");
        System.out.println(Arrays.toString(resultado));
    }
    
    public static double[] jacobi(double[][] coeficientes, double[] soluciones, double tolerancia, int maxIteraciones) {
        int n = soluciones.length;
        double[] x = new double[n];
        double[] xAnterior = new double[n];
        int iteraciones = 0;
        
        while (iteraciones < maxIteraciones) {
            for (int i = 0; i < n; i++) {
                xAnterior[i] = x[i];
            }
            
            for (int i = 0; i < n; i++) {
                double sum = soluciones[i];
                for (int j = 0; j < n; j++) {
                    if (j != i) {
                        sum -= coeficientes[i][j] * xAnterior[j];
                    }
                }
                x[i] = sum / coeficientes[i][i];
            }
            
            double error = 0;
            for (int i = 0; i < n; i++) {
                error = Math.max(error, Math.abs(x[i] - xAnterior[i]));
            }
            
            if (error < tolerancia) {
                break;
            }
            
            iteraciones++;
        }
        
        return x;
    }
}
## Ejercicio 4:
public class JacobiMethod {

    // Función para imprimir una matriz
    public static void printMatrix(double[][] matrix) {
        int n = matrix.length;
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n + 1; j++) {
                System.out.print(matrix[i][j] + " ");
            }
            System.out.println();
        }
    }

    // Función para calcular la norma de la diferencia entre dos vectores
    public static double norm(double[] v1, double[] v2) {
        int n = v1.length;
        double sum = 0;
        for (int i = 0; i < n; i++) {
            sum += Math.pow((v1[i] - v2[i]), 2);
        }
        return Math.sqrt(sum);
    }

    // Función para resolver el sistema de ecuaciones utilizando el método de Jacobi
    public static double[] jacobi(double[][] matrix, int maxIterations, double tolerance) {
        int n = matrix.length;
        double[] x = new double[n];
        double[] x_new = new double[n];
        int iterations = 0;

        while (iterations < maxIterations) {
            for (int i = 0; i < n; i++) {
                double sum = matrix[i][n];
                for (int j = 0; j < n; j++) {
                    if (j != i) {
                        sum -= matrix[i][j] * x[j];
                    }
                }
                x_new[i] = sum / matrix[i][i];
            }

            if (norm(x, x_new) < tolerance) {
                System.arraycopy(x_new, 0, x, 0, n);
                break;
            }

            System.arraycopy(x_new, 0, x, 0, n);
            iterations++;
        }

        return x;
    }

    public static void main(String[] args) {
        double[][] matrix = {
            {12, 3, -5, 1},
            {3, -8, -3, 1},
            {5, 4, -12, 3}
        };

        int maxIterations = 1000;
        double tolerance = 1e-6;

        double[] solution = jacobi(matrix, maxIterations, tolerance);

        System.out.println("Solución:");
        for (int i = 0; i < solution.length; i++) {
            System.out.println("x" + (i + 1) + " = " + solution[i]);
        }
    }
}

## Ejercicio 5:
import java.util.Arrays;

public class MetodoJacobiEjemplo5 {

    public static void main(String[] args) {
        // Definir el nuevo sistema de ecuaciones
        double[][] coeficientes = {
                {5, 2, -3},
                {1, -4, 1},
                {2, 1, -5}
        };

        double[] soluciones = {11, -2, 9}; // Cambiar las soluciones según el nuevo sistema

        // Resolver el sistema utilizando el método de Jacobi
        double[] resultado = jacobi(coeficientes, soluciones, 0.0001, 100);

        // Imprimir las soluciones
        System.out.println("Soluciones:");
        System.out.println(Arrays.toString(resultado));
    }

    public static double[] jacobi(double[][] coeficientes, double[] soluciones, double tolerancia, int maxIteraciones) {
        int n = soluciones.length;
        double[] x = new double[n];
        double[] xAnterior = new double[n];
        int iteraciones = 0;

        while (iteraciones < maxIteraciones) {
            for (int i = 0; i < n; i++) {
                xAnterior[i] = x[i];
            }

            for (int i = 0; i < n; i++) {
                double sum = soluciones[i];
                for (int j = 0; j < n; j++) {
                    if (j != i) {
                        sum -= coeficientes[i][j] * xAnterior[j];
                    }
                }
                x[i] = sum / coeficientes[i][i];
            }

            double error = 0;
            for (int i = 0; i < n; i++) {
                error = Math.max(error, Math.abs(x[i] - xAnterior[i]));
            }

            if (error < tolerancia) {
                break;
            }

            iteraciones++;
        }

        return x;
    }
}

# Método Gauss-Jordan
Este método debe su nombre a Carl Friedrich Gauss y a Wilhelm jordan. Se trata de una serie de algoritmos del algebra lineal para determinar los resultados de un sistema de ecuaciones lineales y así hallar matrices e inversas. El sistema de Gauss se utiliza para resolver un sistema de ecuaciones y obtener las soluciones por medio de la reducción del sistema dado a otro que sea equivalente en el cual cada una de las ecuaciones tendrá una incógnita menos que la anterior. La matriz que resulta de este proceso lleva el nombre que se conoce como forma escalonada.

Este método, permite resolver hasta 20 ecuaciones simultáneas. Lo que lo diferencia del método Gaussiano es que cuando es eliminada una incógnita, se eliminará de todas las ecuaciones restantes, o sea, las que anteceden a la ecuación principal así como de las que la siguen a continuación. De esta manera el paso de eliminación forma una matriz identidad en vez de una matriz triangular. No es necesario entonces utilizar la sustitución hacia atrás para conseguir la solución.
## Algoritmo
Proceso Gauss_Jordan()
    // Definir variables
    Definir matriz[3, 4] como Real
    Definir filas como Entero
    Definir columnas como Entero
    
    // Inicializar la matriz de coeficientes y términos independientes
    matriz[1][1] = 2
    matriz[1][2] = 1
    matriz[1][3] = -1
    matriz[1][4] = 8
    matriz[2][1] = -3
    matriz[2][2] = -1
    matriz[2][3] = 2
    matriz[2][4] = -11
    matriz[3][1] = -2
    matriz[3][2] = 1
    matriz[3][3] = 2
    matriz[3][4] = -3
    
    // Obtener dimensiones de la matriz
    filas = 3
    columnas = 4
    
    // Bucle principal del método de Gauss-Jordan
    Para i desde 1 hasta filas
        // Encontrar el pivote
        pivote = matriz[i][i]
        
        // Dividir la fila por el pivote para hacer el pivote igual a 1
        Para j desde i hasta columnas
            matriz[i][j] = matriz[i][j] / pivote
        FinPara
        
        // Hacer ceros debajo del pivote
        Para k desde 1 hasta filas
            Si (k ≠ i) Entonces
                factor = matriz[k][i]
                Para j d

## Ejercicio 1:
public class GaussJordan {

 public static void main(String[] args) {
        double[][] matriz = {{2, 1, -1, 8},
                              {-3, -1, 2, -11},
                              {-2, 1, 2, -3}};
        
        double[] soluciones = resolverSistema(matriz);
        
        System.out.println("Soluciones:");
        for (int i = 0; i < soluciones.length; i++) {
            System.out.println("x" + (i+1) + " = " + soluciones[i]);
        }
    }
    
    public static double[] resolverSistema(double[][] matriz) {
        int n = matriz.length;
        double[] soluciones = new double[n];
        
        // Proceso de Gauss-Jordan
        for (int i = 0; i < n; i++) {
            // Pivote actual
            double pivote = matriz[i][i];
            
            // Hacer el pivote igual a 1
            for (int j = 0; j < n + 1; j++) {
                matriz[i][j] /= pivote;
            }
            
            // Eliminación hacia adelante y hacia atrás
            for (int k = 0; k < n; k++) {
                if (k != i) {
                    double factor = matriz[k][i];
                    for (int j = 0; j < n + 1; j++) {
                        matriz[k][j] -= factor * matriz[i][j];
                    }
                }
            }
        }
        
        // Extraer las soluciones
        for (int i = 0; i < n; i++) {
            soluciones[i] = matriz[i][n];
        }
        
        return soluciones;
    }
}

## Ejercicios 2:

    public static double[] resolverSistema(double[][] matriz) {
        int n = matriz.length;

        for (int i = 0; i < n; i++) {
            // Hacer que el pivote sea 1
            double pivote = matriz[i][i];
            for (int j = i; j < n + 1; j++) {
                matriz[i][j] /= pivote;
            }

            // Hacer 0 en los otros elementos de la columna del pivote
            for (int k = 0; k < n; k++) {
                if (k != i) {
                    double factor = matriz[k][i];
                    for (int j = i; j < n + 1; j++) {
                        matriz[k][j] -= factor * matriz[i][j];
                    }
                }
            }
        }

        // Extraer la soluciÃ³n
        double[] solucion = new double[n];
        for (int i = 0; i < n; i++) {
            solucion[i] = matriz[i][n];
        }

        return solucion;
    }
## Ejercicios 3:
public class GaussJordan {

    public static void gaussJordan(double[][] matriz) {
        int filas = matriz.length;
        int columnas = matriz[0].length;

        // Bucle principal del método de Gauss-Jordan
        for (int i = 0; i < filas; i++) {
            // Encontrar el pivote
            double pivote = matriz[i][i];

            // Dividir la fila por el pivote para hacer el pivote igual a 1
            for (int j = i; j < columnas; j++) {
                matriz[i][j] /= pivote;
            }

            // Hacer ceros debajo del pivote
            for (int k = 0; k < filas; k++) {
                if (k != i) {
                    double factor = matriz[k][i];
                    for (int j = i; j < columnas; j++) {
                        matriz[k][j] -= factor * matriz[i][j];
                    }
                }
            }
        }

        // Imprimir las soluciones
        System.out.println("Soluciones:");
        for (int i = 0; i < filas; i++) {
            System.out.println("x" + (i + 1) + ": " + matriz[i][columnas - 1]);
        }
    }

    public static void main(String[] args) {
        // Definir la matriz de coeficientes y términos independientes
        double[][] matriz = {
                {2, -6, 4, 6},
                {-3, -1, -3, 1},
                {-4, -1, 0, 6}
        };

        // Resolver el sistema de ecuaciones lineales
        gaussJordan(matriz);
    }
}
## Ejercicio 4:
public class GaussJordanEjemplo4 {

    // Ecuacion a Resolver
    // 4x - 2y +  z =  9
    // 3x -  y + 4z = 10
    // 2x +  y -  z =  8

    public static void main(String[] args) {
        double[][] matriz = {
                {4, -2, 1, 9},
                {3, -1, 4, 10},
                {2, 1, -1, 8}
        };

        resolverSistema(matriz);
    }

    public static void resolverSistema(double[][] matriz) {
        int filas = matriz.length;
        int columnas = matriz[0].length;

        // Escalonar la matriz
        for (int i = 0; i < filas; i++) {
            double pivote = matriz[i][i];
            for (int j = i + 1; j < filas; j++) {
                double factor = matriz[j][i] / pivote;
                for (int k = i; k < columnas; k++) {
                    matriz[j][k] -= matriz[i][k] * factor;
                }
            }
        }

        // Reducir la matriz a forma escalonada reducida
        for (int i = filas - 1; i >= 0; i--) {
            double pivote = matriz[i][i];
            for (int j = i - 1; j >= 0; j--) {
                double factor = matriz[j][i] / pivote;
                for (int k = i; k < columnas; k++) {
                    matriz[j][k] -= matriz[i][k] * factor;
                }
            }
        }

        // Normalizar la matriz
        for (int i = 0; i < filas; i++) {
            double pivote = matriz[i][i];
            for (int j = i; j < columnas; j++) {
                matriz[i][j] /= pivote;
            }
        }

        // Imprimir la solución
        for (int i = 0; i < filas; i++) {
            System.out.print("x" + (i + 1) + " = " + matriz[i][columnas - 1]);
            System.out.println();
        }
    }
}
## Ejercicio 5:
public class GaussJordanEjemplo5 {

    // Ecuacion a Resolver
    // 3x -  y + 4z = 7
    // 2x +  y - 3z = 8
    //  x - 4y + 2z = 9

    public static void main(String[] args) {
        double[][] matriz = {
                {3, -1, 4, 7},
                {2, 1, -3, 8},
                {1, -4, 2, 9}
        };

        resolverSistema(matriz);
    }

    public static void resolverSistema(double[][] matriz) {
        int filas = matriz.length;
        int columnas = matriz[0].length;

        // Escalonar la matriz
        for (int i = 0; i < filas; i++) {
            double pivote = matriz[i][i];
            for (int j = i + 1; j < filas; j++) {
                double factor = matriz[j][i] / pivote;
                for (int k = i; k < columnas; k++) {
                    matriz[j][k] -= matriz[i][k] * factor;
                }
            }
        }

        // Reducir la matriz a forma escalonada reducida
        for (int i = filas - 1; i >= 0; i--) {
            double pivote = matriz[i][i];
            for (int j = i - 1; j >= 0; j--) {
                double factor = matriz[j][i] / pivote;
                for (int k = i; k < columnas; k++) {
                    matriz[j][k] -= matriz[i][k] * factor;
                }
            }
        }

        // Normalizar la matriz
        for (int i = 0; i < filas; i++) {
            double pivote = matriz[i][i];
            for (int j = i; j < columnas; j++) {
                matriz[i][j] /= pivote;
            }
        }

        // Imprimir la solución
        for (int i = 0; i < filas; i++) {
            System.out.print("x" + (i + 1) + " = " + matriz[i][columnas - 1]);
            System.out.println();
        }
    }
}

# Eliminación Gaussiana
La eliminación gaussiana es el proceso de usar operaciones de fila válidas en una matriz hasta que esté en forma escalonada de fila reducida. Hay tres tipos de operaciones de fila válidas que se pueden realizar en una matriz.
OP1: intercambia dos filas.
OP2: multiplica todas las entradas de una fila por un número distinto de cero.
OP3: agregue un múltiplo de una fila a una fila de destino. (Nota: la fila de destino es la única fila que se cambia en este proceso).
Es importante darse cuenta de que estas son solo las reglas del juego. La forma en que los apliquemos en una situación determinada dependerá de la matriz que se nos dé. Tenga en cuenta que nuestro objetivo es transformar la matriz en una forma más simple, llamada forma escalonada de fila reducida (RREF) , utilizando una serie de estas operaciones de fila.
## Algoritmo
Proceso Eliminacion_Gaussiana()
    // Definir variables
    Definir matriz[3, 4] como Real
    Definir n como Entero
    Definir i, j, k como Entero
    Definir pivote, factor como Real

    // Inicializar la matriz de coeficientes
    matriz[1][1] = 2
    matriz[1][2] = 1
    matriz[1][3] = -1
    matriz[1][4] = 8
    matriz[2][1] = -3
    matriz[2][2] = -1
    matriz[2][3] = 2
    matriz[2][4] = -11
    matriz[3][1] = -2
    matriz[3][2] = 1
    matriz[3][3] = 2
    matriz[3][4] = -3

    // Obtener el tamaño de la matriz
    n = 3

    // Aplicar eliminación gaussiana
    Para k desde 1 hasta n-1
        Para i desde k+1 hasta n
            factor = matriz[i][k] / matriz[k][k]
            Para j desde k+1 hasta n+1
                matriz[i][j] = matriz[i][j] - factor * matriz[k][j]
            FinPara
        FinPara
    FinPara

    // Resolver el sistema triangular superior
    Para i desde n hasta 1
        pivote = matriz[i][i]
        Para j desde n+1 hasta n+1
            matriz[i][j] = matriz[i][j] / pivote
        FinPara
        Para k desde i-1 hasta 1
            factor = matriz[k][i]
            Para j desde n+1 hasta n+1
                matriz[k][j] = matriz[k][j] - factor * matriz[i][j]
            FinPara
        FinPara
    FinPara

    // Imprimir las soluciones
    Escribir("Soluciones:")
    Para i desde 1 hasta n
        Escribir("x" + i + ": " + matriz[i][n+1])
    FinPara
FinProceso

## Ejercicio 1:
package EliminacionGaussiana;

public class EliminacionGaussiana {

    public static void main(String[] args) {
        double[][] sistemaEcuaciones = {
                {2, 1, -1, 8},
                {-3, -1, 2, -11},
                {-2, 1, 2, -3}
        };

        resolverSistema(sistemaEcuaciones);
    }

    public static void resolverSistema(double[][] matriz) {
        int filas = matriz.length;
        int columnas = matriz[0].length - 1; // Ignoramos la última columna (soluciones)

        // Aplicar eliminación gaussiana
        for (int i = 0; i < filas - 1; i++) {
            for (int j = i + 1; j < filas; j++) {
                double factor = matriz[j][i] / matriz[i][i];
                for (int k = i; k < columnas + 1; k++) {
                    matriz[j][k] -= factor * matriz[i][k];
                }
            }
        }

        // Realizar sustitución hacia atrás
        double[] soluciones = new double[filas];
        for (int i = filas - 1; i >= 0; i--) {
            double suma = 0;
            for (int j = i + 1; j < filas; j++) {
                suma += matriz[i][j] * soluciones[j];
            }
            soluciones[i] = (matriz[i][columnas] - suma) / matriz[i][i];
        }

        // Mostrar las soluciones
        System.out.println("Soluciones:");
        for (int i = 0; i < filas; i++) {
            System.out.println("x" + (i + 1) + " = " + soluciones[i]);
        }
    }
}
## Ejercicio 2:

    public static void main(String[] args) {
        double[][] matriz = {{2, 1, -1, 8},
                              {-3, -1, 2, -11},
                              {-2, 1, 2, -3}};
        
        double[] soluciones = resolverSistema(matriz);
        
        System.out.println("Soluciones:");
        for (int i = 0; i < soluciones.length; i++) {
            System.out.println("x" + (i+1) + " = " + soluciones[i]);
        }
    }
    
    public static double[] resolverSistema(double[][] matriz) {
        int n = matriz.length;
        double[] soluciones = new double[n];
        
        // Eliminación hacia adelante
        for (int i = 0; i < n; i++) {
            for (int k = i+1; k < n; k++) {
                double factor = matriz[k][i] / matriz[i][i];
                for (int j = i; j < n+1; j++) {
                    matriz[k][j] -= factor * matriz[i][j];
                }
            }
        }
        
        // Sustitución hacia atrás
        for (int i = n-1; i >= 0; i--) {
            soluciones[i] = matriz[i][n] / matriz[i][i];
            for (int k = i-1; k >= 0; k--) {
                matriz[k][n] -= matriz[k][i] * soluciones[i];
            }
        }
        
        return soluciones;
    }
}

## Ejercicio 3:
package eliminaciongaussiana;

import java.util.Scanner;
import javax.swing.JOptionPane;

public class EliminacionGaussiana {

    public static void main(String[] args) {
JOptionPane.showMessageDialog(null, "Ecuaciones por Eliminación Gaussiana");

        int n = Integer.parseInt(JOptionPane.showInputDialog("Ingrese el número de ecuaciones: "));
        float[][] matriz = new float[n][n + 1];

        JOptionPane.showMessageDialog(null, "Ingrese los coeficientes y resultados de cada ecuación:");

        for (int i = 0; i < n; i++) {
            String[] coeficientes = JOptionPane.showInputDialog("Ecuación " + (i + 1) + ":").split(" ");
            for (int j = 0; j < n + 1; j++) {
                matriz[i][j] = Float.parseFloat(coeficientes[j]);
            }
        }

        for (int i = 0; i < n - 1; i++) {
            for (int k = i + 1; k < n; k++) {
                float factor = matriz[k][i] / matriz[i][i];
                for (int j = i; j < n + 1; j++) {
                    matriz[k][j] -= factor * matriz[i][j];
                }
            }
        }

        float[] soluciones = new float[n];
        for (int i = n - 1; i >= 0; i--) {
            float suma = 0;
            for (int j = i + 1; j < n; j++) {
                suma += matriz[i][j] * soluciones[j];
            }
            soluciones[i] = (matriz[i][n] - suma) / matriz[i][i];
        }

        StringBuilder mensaje = new StringBuilder("Soluciones:\n");
        for (int i = 0; i < n; i++) {
            mensaje.append("x[").append(i + 1).append("] = ").append(soluciones[i]).append("\n");
        }
        JOptionPane.showMessageDialog(null, mensaje.toString());
    }
}
## Ekercicio 4:
public class EliminacionGaussianaEjemplo5 {

    public static void main(String[] args) {
        double[][] matriz = {{3, -1, 4, 7},
                {2, 1, -3, 8},
                {1, -4, 2, 9}};

        double[] soluciones = resolverSistema(matriz);

        System.out.println("Soluciones:");
        for (int i = 0; i < soluciones.length; i++) {
            System.out.println("x" + (i+1) + " = " + soluciones[i]);
        }
    }

    public static double[] resolverSistema(double[][] matriz) {
        int n = matriz.length;
        double[] soluciones = new double[n];

        // Eliminación hacia adelante
        for (int i = 0; i < n; i++) {
            for (int k = i+1; k < n; k++) {
                double factor = matriz[k][i] / matriz[i][i];
                for (int j = i; j < n+1; j++) {
                    matriz[k][j] -= factor * matriz[i][j];
                }
            }
        }

        // Sustitución hacia atrás
        for (int i = n-1; i >= 0; i--) {
            soluciones[i] = matriz[i][n] / matriz[i][i];
            for (int k = i-1; k >= 0; k--) {
                matriz[k][n] -= matriz[k][i] * soluciones[i];
            }
        }

        return soluciones;
    }
}
## Ejercicio 5:
public class EliminacionGaussianaEjemplo4 {

    public static void main(String[] args) {
        double[][] matriz = {{4, -2, 1, 9},
                {3, -1, 4, 10},
                {2, 1, -1, 8}};

        double[] soluciones = resolverSistema(matriz);

        System.out.println("Soluciones:");
        for (int i = 0; i < soluciones.length; i++) {
            System.out.println("x" + (i+1) + " = " + soluciones[i]);
        }
    }

    public static double[] resolverSistema(double[][] matriz) {
        int n = matriz.length;
        double[] soluciones = new double[n];

        // Eliminación hacia adelante
        for (int i = 0; i < n; i++) {
            for (int k = i+1; k < n; k++) {
                double factor = matriz[k][i] / matriz[i][i];
                for (int j = i; j < n+1; j++) {
                    matriz[k][j] -= factor * matriz[i][j];
                }
            }
        }

        // Sustitución hacia atrás
        for (int i = n-1; i >= 0; i--) {
            soluciones[i] = matriz[i][n] / matriz[i][i];
            for (int k = i-1; k >= 0; k--) {
                matriz[k][n] -= matriz[k][i] * soluciones[i];
            }
        }

        return soluciones;
    }
}

# INTEGRANTES: Equipo 2
López Gutierrez Nili Estefanía
Mogollan Agis Kate Michelle
Hernández Zavala Gabriel
De Santiago Castillo Evelyn
Alvarez Rivera Angel
Ruiz Aguilera Diego 
