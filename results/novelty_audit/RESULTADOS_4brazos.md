# MFIOP — Prueba de 4 brazos sobre 5SYF y 3J6F
## Control de novedad — PRELIMINAR, NO verificado  — 2026-10-09

## ESTADO: resultado sugestivo y reproducible, pero NO verificado.
No usar como "demostrado" sin completar las 3 condiciones de verificacion (abajo).

## Lo que SÍ sostiene el artefacto (reproducible)
| Brazo | Numero | Que dice |
|-------|--------|----------|
| 1. MFIOP puro | 5SYF z=+1.19 (pct 85); 3J6F z=+0.92 (pct 81) vs nulo orient | Hay senal sobre azar |
| 2. Topologico | ratio sin orient 1.478 vs con orient 1.639 | Orientacion aporta ~11% |
| 3. amida-I (simplificado) | amida-I ratio 1.027 ; MFIOP(Trp) ratio 1.639 | Indicio de eje independiente |
| 4. Nulo geometrico | ambas superan el nulo de orientacion | Estructura real > azar |

## LO QUE NO ESTA VERIFICADO (por honestidad - no afirmar lo contrario)
1. El brazo 3 compara UN ratio por sistema, SIN barras de error ni perturbacion.
   No es un test estadistico; es una comparacion de un punto.
2. El "amida-I" usado es SIMPLIFICADO: mismo kappa/r^3 con dipolos unitarios sobre C=O,
   NO la parametrizacion real de amida-I (mapas Jansen/Torii-Tasumi, energias de sitio,
   acoplamiento nearest-neighbor por Ramachandran). Es "MFIOP sobre C=O", no el amida-I
   de la literatura. Un referee de espectroscopia exigiria el modelo real.
3. "range" (max-min autovalor) es sensible a outliers y al tamano N (72 vs 8 Trp).
   Falta metrica normalizada por N (p.ej. participation ratio) para comparacion justa.

## Para pasar de "sugestivo" a "verificado" faltan 3 cosas
(a) amida-I con parametrizacion real de la literatura.
(b) barras de error via perturbacion geometrica/orientacional en el diseno comparativo.
(c) metrica normalizada por N.

## Formulacion HONESTA (lo maximo defendible HOY)
"Un analisis preliminar sugiere que el descriptor de triptofanos captura variacion
conformacional entre arquitecturas de tubulina que un modelo vibracional simplificado
del backbone no refleja; la verificacion requiere un control amida-I completo y
analisis de error." -- NO afirmar "demostrado independiente" todavia.

## Scripts: descriptor.py (brazos 1+4), topo.py (brazo 2), amide3.py (brazo 3 simplificado)
