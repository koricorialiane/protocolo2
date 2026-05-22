# Protocolo Titán · versión profesional

Esta versión transforma la app en un panel técnico de ingeniería para justificar una red GSM/EDGE táctica.

## Ejecutar

```bash
pip install -r requirements.txt
streamlit run src/protocolo_titan/streamlit_app.py
```

## Mejoras principales

- Dashboard Streamlit profesional con sidebar de parámetros, pestañas y exportación CSV.
- Escenario A:
  - curva continua Doppler-velocidad;
  - curva continua de tiempo de coherencia frente al timeslot GSM;
  - rotación de fase por timeslot;
  - traza de fading tipo Jakes para Rayleigh y Rician;
  - presupuesto de enlace con FSPL y Okumura-Hata urbano;
  - margen de enlace en borde de celda.
- Escenario B:
  - asignación ARFCN uplink/downlink GSM-900;
  - BCCH/TCH con política de hopping y potencia;
  - distancia de reutilización `D = R sqrt(3N)`;
  - SIR co-canal de primer anillo;
  - capacidad TCH aproximada;
  - probabilidad de bloqueo Erlang B;
  - curva de ruido integrada según RBW;
  - compromiso sensibilidad-tiempo de barrido;
  - checklist RED.

## Archivos modificados

- `src/protocolo_titan/streamlit_app.py`
- `src/protocolo_titan/ui_charts.py`
- `src/protocolo_titan/propagation.py`
- `src/protocolo_titan/scenario_a.py`
- `src/protocolo_titan/cellular_planning.py`
- `src/protocolo_titan/scenario_b.py`
- `src/protocolo_titan/instrumentation.py`

## Nota técnica

Las curvas de fading se han hecho con una suma de osciladores tipo Jakes para que el canal varíe de forma correlada y dependiente del Doppler. No sustituye a un simulador radio completo, pero es mucho más defendible que una curva plana o una aleatoriedad sin física.
