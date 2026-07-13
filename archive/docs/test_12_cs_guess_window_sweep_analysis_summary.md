# Test 12: análisis del barrido de `cs_guess` y ventana efectiva

Este documento resume el análisis preparado para Test 12. Los valores numéricos finales deben llenarse ejecutando:

```matlab
experiments/analysis/analyze_test_12_level01_model_comparison.m
experiments/analysis/analyze_test_12_level02_grouped_generalization.m
experiments/analysis/analyze_test_12_level03_compare_test11.m
```

## 1. Qué se hizo en Test 12

Test 12 extiende Test 11 variando explícitamente `REQ_cs_guess`, además de `REQ_M`, `SIM_cs_bg`, `SIM_f0` y el modelo de onda. El objetivo es desacoplar el tamaño nominal de ventana de REQ del tamaño efectivo de ventana relativo a la longitud de onda verdadera.

## 2. Por qué se varió `REQ_cs_guess`

En Test 11, `cs_guess` estaba fijo. Eso acoplaba artificialmente:

```matlab
REQ_M
REQ_cs_guess
SIM_cs_bg
M_eff_true_diag
```

Test 12 permite preguntar si el modelo falla por `REQ_M` nominal, por el valor usado como `cs_guess`, o por el tamaño efectivo real de la ventana.

## 3. Modelos evaluados

Los modelos operacionales son:

- `LocalOnly`: usa solo features locales de cada ventana.
- `GlobalOnly`: usa solo features globales del campo completo.
- `HybridLocalGlobal`: combina features locales y globales.

El análisis también permite un modelo diagnóstico con `M_eff_true_diag`, pero este no debe usarse como conclusión operacional porque depende de la velocidad verdadera simulada.

## 4. Feature sets comparados

Los feature sets son:

- `NoCsGuess`: equivalente al enfoque de Test 11, sin `REQ_cs_guess`, `M_eff_guess` ni `M_eff_true_diag`.
- `WithCsGuess`: agrega `REQ_cs_guess`.
- `WithMeffGuess`: agrega `M_eff_guess`.
- `WithCsGuessAndMeffGuess`: agrega `REQ_cs_guess` y `M_eff_guess`.
- `DiagnosticWithMeffTrue`: agrega `M_eff_true_diag`; solo diagnóstico.

## 5. ¿`REQ_cs_guess` mejora?

TODO: ejecutar Level 01 y Level 02. La respuesta debe salir de:

```text
outputs/test_12_cs_guess_window_sweep/.../analysis/level_01_model_comparison/tables/level12_level01_delta_mape_vs_NoCsGuess.csv
outputs/test_12_cs_guess_window_sweep/.../analysis/level_02_grouped_generalization/tables/level12_level02_delta_mape_vs_NoCsGuess.csv
```

Delta negativo significa mejora frente a `NoCsGuess`.

## 6. ¿`M_eff_guess` mejora?

TODO: verificar en las mismas tablas de delta. La comparación principal debe mirar:

- `WithMeffGuess - NoCsGuess`
- `WithCsGuessAndMeffGuess - NoCsGuess`

## 7. ¿Test 12 mejora respecto a Test 11?

TODO: ejecutar Level 03. La comparación global queda en:

```text
level12_level03_test11_vs_test12_global_metrics.csv
```

La comparación matched/subset queda en:

```text
level12_level03_test11_vs_test12_matched_metrics.csv
```

## 8. Comparación global vs matched

La comparación global usa las métricas completas de Test 11 y Test 12. Sirve como referencia general, pero no es completamente apples-to-apples porque Test 12 cambió el diseño experimental.

La comparación matched restringe los datos a condiciones comunes cuando es posible:

- valores comunes de `SIM_f0`;
- valores comunes de `SIM_cs_bg`;
- valores comunes de `REQ_M`;
- valores comunes de `SIM_WaveModel`;
- `REQ_cs_guess == 3.0` en Test 12, si existe.

Esta segunda comparación es la más importante para decidir si Test 12 mejora de forma justa frente a Test 11.

## 9. Limitaciones pendientes

- Test 12 todavía no analiza SNR.
- Test 12 todavía no incluye atenuación.
- Test 12 no reemplaza validaciones heterogéneas como bilayer o k-Wave.
- `M_eff_true_diag` es una variable oracle y solo debe usarse para diagnóstico.
- La utilidad de `M_eff_guess` debe evaluarse con grouped splits, especialmente leave-one-`REQ_cs_guess` y leave-one-`SIM_cs_bg`.

## 10. Próximos pasos sugeridos

1. Correr Level 01 para cerrar la ablation dentro de Test 12.
2. Correr Level 02 para grouped generalization con `REQ_cs_guess` y `SIM_cs_bg` no vistos.
3. Correr Level 03 para comparar Test 12 contra Test 11 global y matched.
4. Crear Test 13 para SNR.
5. Crear otro test para SNR dependiente de profundidad.
6. Agregar atenuación.
7. Volver a bilayer.
8. Validar después con k-Wave.
