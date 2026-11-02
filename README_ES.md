# DE_Particle_paper — guía breve en español

Proyecto reproducible para reconstrucción inversa de un sector escalar 2HDM + singlete complejo mediante Differential Evolution. Contiene **3000 optimizaciones nuevas**, 1000 tripletes objetivo y tres reinicios por target (semillas base 11, 23 y 37), con `v=246.22 GeV`.

El corte final es `|z_i|<1` para **cada uno de los cinco residuos**: tres masas, masa escalar cargada y acoplamiento κ. Se aceptan **2971/3000** ejecuciones y los 1000 targets tienen una o más soluciones aceptadas. Las σ son escalas numéricas adoptadas en la reconstrucción, no necesariamente errores experimentales.

## Uso

Abre la carpeta raíz del repositorio en VS Code, crea un entorno con Python 3.12 e instala `requirements.txt`:

```powershell
python -m venv .venv
.venv\Scripts\Activate.ps1
python -m pip install -r requirements.txt
python scripts/verify_dataset.py
```

Abre `notebooks/00_start_here.ipynb` y ejecuta sus celdas. Los seis notebooks están escritos en inglés para facilitar su utilización internacional, llevan instrucciones en Markdown, usan rutas relativas y conservan salidas pequeñas para visualizarlas directamente en GitHub. No necesitas repetir las 3000 optimizaciones para ejecutar el análisis.

Para una reproducción exhaustiva del entorno original, consulta `data/raw/production_environment_full.txt` (inventario completo del equipo Windows, convertido a UTF-8 para visualizarlo en GitHub, **no** lista mínima recomendada para instalar). La versión recomendada de dependencias es `requirements.txt`.

**Publicación:** antes de habilitar visibilidad pública, los tres autores deben aprobar los datos, la licencia y las citas; véase `docs/publication_checklist.md`. El conector de GitHub vinculado permite que más adelante subamos los archivos cuando exista el repositorio `DE_Particle_paper`.
