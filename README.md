## Bem vindo ao repositório do código de localização cilindrica

Este código pode ser acessado no [GitHub](https://github.com/LeonardoSirino/CylindricalLocation).

### Descrição

Este repositório tem por objetivo fornecer uma biblioteca para a localização de fontes de Emissão Acústica em vasos de pressão com tampo elispoidal. Para isso são utilizados métodos de minimização e inteligência artificial que visam melhorar a precisão das técnica atualmente utilizadas.

Outro diferencial é a utilização do cálculo de geodésicas em elipsoides, desta forma, não existem distorções nas distâncias, algo comum nas técnicas tradicionais.

```markdown
from Routines.CilLoc import *

Locate = CylindricalLocation(1000.0, 1000.0, 750 * m.pi / 4)
CoordList = [{"Xcord": 200.0, "Ycord": 1200.0},
             {"Xcord": 100.0, "Ycord": 900.0},
             {"Xcord": 200.0, "Ycord": 700.0},
             {"Xcord": 2500.0, "Ycord": -300}]

for coord in CoordList:
    Locate.AddSensor(Xcord=coord.get("Xcord"), Ycord=coord.get("Ycord"))

Locate.calcAllDist(600, 900, True)
```
