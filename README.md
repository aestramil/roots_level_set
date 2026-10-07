# roots_level_set
# Códigos de la tesis de Maestría

Este repositorio contiene los códigos utilizados en los experimentos
numéricos del Capítulo 3 de la tesis

**Dos aplicaciones de las cópulas: búsqueda estocástica de conjuntos de
nivel y análisis de datos educativos**

correspondiente a la Maestría en Ingeniería Matemática de la
Universidad de la República.

**Autor:** Agustín Estramil  
**Directores:** Leonardo Moreno y Emilien Joly  
**Año:** 2026

## Descripción

El Capítulo 3 estudia un procedimiento estocástico para la búsqueda de
raíces y la aproximación de conjuntos de nivel basado en propuestas
guiadas por cópulas.

Los códigos incluidos en este repositorio permiten reproducir los
principales experimentos numéricos presentados en dicho capítulo. El
material se encuentra organizado de acuerdo con las distintas etapas
del estudio y con la aplicación a datos de *lifetime* de lingotes de
silicio.

## Estructura del repositorio

### `01_caso_univariado/`

Códigos correspondientes a los experimentos de búsqueda de raíces en
una dimensión.

Esta carpeta contiene las implementaciones y funciones de prueba
utilizadas para estudiar el comportamiento del procedimiento basado en
cópulas en el caso univariado.

### `02_conjuntos_nivel/`

Códigos correspondientes a la extensión del procedimiento a la
aproximación de conjuntos de nivel en dimensión dos.

Incluye las implementaciones utilizadas para la generación de
candidatos mediante cópulas y el mecanismo de Metropolis--Hastings,
así como los procedimientos de refinamiento considerados en la tesis.

### `03_conjuntos_nivel_replicas/`

Códigos utilizados para realizar las réplicas de los experimentos
bidimensionales y obtener las medidas de desempeño reportadas en la
tesis.

Esta carpeta contiene las implementaciones empleadas para evaluar el
comportamiento de los distintos métodos sobre múltiples realizaciones
de los experimentos.

### `04_ejemplo_datos_silicio/`

Códigos correspondientes a la aplicación del procedimiento a datos de
*lifetime* en lingotes de silicio.

El objetivo de esta aplicación es aproximar el conjunto de nivel
considerado en la tesis a partir de observaciones disponibles sobre
una grilla espacial.

## Datos de la aplicación a silicio

Los archivos

- `data3`
- `data4`
- `github_datos_silicio`

contienen información utilizada en la aplicación a los datos de
*lifetime* presentada en el Capítulo 3.

Los scripts correspondientes al análisis de estos datos se encuentran
en la carpeta `04_ejemplo_datos_silicio/`.

## Software

Los experimentos fueron implementados en **R**.

Los paquetes requeridos se especifican en los scripts correspondientes
a cada experimento. Estos incluyen herramientas para simulación y
modelización mediante cópulas, procesos gaussianos, cálculo numérico
y visualización gráfica.

## Reproducción de los experimentos

La organización del repositorio sigue el orden general de los
experimentos presentados en el Capítulo 3 de la tesis.

Los scripts de cada carpeta pueden ejecutarse siguiendo las
indicaciones incluidas en el código. Cuando corresponde, las semillas
utilizadas para la generación pseudoaleatoria se encuentran fijadas
con el objetivo de facilitar la reproducibilidad de los resultados.

La carpeta `03_conjuntos_nivel_replicas/` contiene específicamente los
códigos empleados para repetir los experimentos bidimensionales y
calcular las medidas de desempeño utilizadas en las comparaciones.

## Correspondencia con la tesis

Para la descripción matemática de los métodos, la definición de las
métricas de desempeño, el diseño de los experimentos y la interpretación
de los resultados, se remite al Capítulo 3 de la tesis.

El repositorio tiene como objetivo complementar dicha descripción
proporcionando los códigos utilizados para generar los resultados
numéricos.

## Contacto

Agustín Estramil  
agustin.estramil@fcea.edu.uy

Instituto de Estadística  
Facultad de Ciencias Económicas y de Administración  
Universidad de la República  
Uruguay
