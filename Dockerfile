FROM python:3.13.9-slim

ENV PYTHONDONTWRITEBYTECODE=1 \
    PYTHONUNBUFFERED=1 \
    MPLBACKEND=Agg \
    PORT=10000

WORKDIR /app

COPY requirements.txt /app/requirements.txt
RUN pip install --no-cache-dir -r /app/requirements.txt

# The web application and only the scientific source files required at runtime.
COPY tempo-design-explorer /app/tempo-design-explorer
COPY oscillator/code/Mechanistic_ODE_Global_Sensitivity_Analysis.py /app/oscillator/code/Mechanistic_ODE_Global_Sensitivity_Analysis.py
COPY oscillator/code/Shared_PLtetO1_Period_Knob_Design_Map.py /app/oscillator/code/Shared_PLtetO1_Period_Knob_Design_Map.py
COPY shutdown_model/shutdown_core.py /app/shutdown_model/shutdown_core.py
COPY rdfmodel_new/model/zhao_core.py /app/rdfmodel_new/model/zhao_core.py

RUN useradd --create-home --uid 10001 tempo
USER tempo

EXPOSE 10000

CMD ["sh", "-c", "python tempo-design-explorer/server.py --host 0.0.0.0 --port ${PORT:-10000}"]
