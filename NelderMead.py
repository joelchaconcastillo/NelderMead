import numpy as n
import time
import sys
import time
import random
# Funcion de evaluacion
def FuncionObjetivo(x):
	z = 0
	for i in range(len(x)-1):
		#z+= x[i]**2
		z+= 100*((x[i+1]-(x[i]**2))**2) + (1-x[i])**2
	return -z
## Esta funcion evalua los puntos del simplejo y regresa
## un vector con los indices ordenados correspondientes a los puntos
def OrdenarPuntos(X, F):
	IndexEval = n.argsort(F)[::-1]
	return IndexEval
##Realiza el calculo del centroide
def PromedioReflejado(X, IndexF):
	Temp = X[IndexF[0:IndexF.size-1],:]
	return n.mean(Temp , axis =0)
##Realiza la reflexion del vertice
def Reflexion(X,Promedio, IndexF, Alfa):
	return Promedio + Alfa*( Promedio - X[IndexF[IndexF.size-1],:])
##Genera un nuevo simplejo dada la informacion con la que se cuenta
def NuevoSimplejo(X, tao, Xp, IndexF):
	for i in range(Xp.size+1):
		if i == 0:
			X[i,:] = X[IndexF[0],:]
		elif i < Xp.size-1:
			X[i,:] =tao*X[ IndexF[0],:] + (1.0-tao)*X[IndexF[i],:]
		else:
			X[i,:] =tao*X[ IndexF[0],:] + (1.0-tao)*Xp
	return X
##Realiza la evaluacion de cada punto (fila) que corresponde al simplejo
##regresa un vector de evaluaciones
def Evaluar(X, Dim):
	F=[]
	for i in range(Dim+1):
		F.append( FuncionObjetivo(X[i,:]) )
	return F
##Es el criterio de paro el cual se establece como la norma
##de la diferencia entre el promedio y el peor punto
def CriterioParo(Promedio, Peor):
	Diferencia = Promedio-Peor
	return n.linalg.norm(Promedio-Peor)
##Es el algoritmo de NelderMead
def NelderMead(X0, Dim, Tolerancia, Maxite):
	##Parametros de configuracion
	Delta=0.5
	Alfa = 1
	Gamma = 1
	Beta = 0.5
	tao = 0.5
	#Se construye la matriz incial copiando el punto inicial en todos los puntos
	#del simplejo y se agrega un ruido a la diagonal de la matriz compuesta por los uevos puntos
	#esto se considera en la matriz temporal
	X = n.eye(Dim, dtype=float)
	Temporal = n.zeros((Dim,Dim),float)
	n.fill_diagonal(Temporal, Delta)
	X[:,:] = X0
	X +=Temporal
	X = n.vstack([X0, X])
	cont = 0
	Promedio = X[0,:]
	IndexF = range(Dim)

	while CriterioParo(Promedio, X[IndexF[Dim-1],:]) > Tolerancia and cont < Maxite:
		F = Evaluar(X, Dim)
		IndexF = OrdenarPuntos(X, F)
		Promedio = PromedioReflejado(X, IndexF)
		XReflexion = Reflexion(X,Promedio, IndexF, Alfa)
		if FuncionObjetivo(XReflexion) > F[IndexF[0]]:
			Expansion = XReflexion + Gamma*(XReflexion - Promedio)
			if FuncionObjetivo(Expansion) > F[IndexF[0]] :
				X[IndexF[Dim], :] = Expansion
			else:
				X[IndexF[Dim], :] = XReflexion
		elif FuncionObjetivo(XReflexion) > F[ IndexF[Dim-1]]:
			X[IndexF[Dim], : ] = XReflexion
		else:
			if FuncionObjetivo(XReflexion) > F[IndexF[Dim]]:
				Xp = XReflexion
			else:
				Xp = X[IndexF[Dim],:]
			Contraccion = Promedio + Beta*(Xp - Promedio)
			if FuncionObjetivo(Contraccion) > F[IndexF[Dim]]:
				X[IndexF[Dim],:] = Contraccion
			else:
				X = NuevoSimplejo(X, tao, Xp, IndexF)
		cont+=1
	return (X[Dim,:], cont)

Dimension = 2
Tolerancia = 1e-9
IteracionesMaximas = 1e86
##Rangos de la distribucion U~(Min, Max)
##que generaran los puntos iniciales de simplejo
Min = -10
Max = 10
X0 = [Min,Max]
##Obtener el promedio de las evaluacion
SumaIteraciones = 0;
SumaEvaluaciones = 0
SumaTiempos = 0
SizeMuestra = 30
for i in range(SizeMuestra):
	t0 = time.clock()
	RMin =random.uniform(Min,Max)
	RMax = random.uniform(Min,Max)
	Resultado = NelderMead(X0, Dimension, Tolerancia,IteracionesMaximas)
	SumaIteraciones += Resultado[1]
	SumaEvaluaciones += FuncionObjetivo(Resultado[0])
	SumaTiempos += time.clock()
print "Promedio iteraciones: ",SumaIteraciones/SizeMuestra
print "Promedio evaluaciones: ",SumaEvaluaciones/SizeMuestra
print "Promedio tiempos: ",SumaTiempos/SizeMuestra
