# huxt-forecast

Kokeellinen CME-saapumisennuste Maahan [HUXt](https://pypi.org/project/huxt/)-aurinkotuulimallilla
(University of Reading). Tulos päivittyy tiedostoon **`docs/huxt_forecast.json`**, jota revontulisivu lukee.

## Miten se toimii

1. **CME:t**: NASA DONKI:n WSA-Enlil-simulaatioiden kartio-CME-syötteet (`cmeInputs`):
   suunta, leveys, nopeus ja aika 21,5 auringonsäteen kohdalla. Mukaan otetaan viimeisen 5 vrk CME:t,
   jotka ovat Maan puolella (|pituus| ≤ 90°) ja voivat vielä olla matkalla.
2. **Taustatuuli**: NOAA SWPC:n mitattu aurinkotuuli (`json/rtsw/rtsw_wind_1m.json`) ja WSA-Enlil-ennuste Maan kohdalla
   (`enlil_time_series.json`). Aikasarja muunnetaan
   Carringtonin pituusasteprofiiliksi ja kuvataan HUXtin sisäreunalle (21,5 Rs).
3. **Ensemble**: oletuksena 500 ajoa. Jokaisessa CME:n suuntaa (±15°), leveyttä (±20 %), nopeutta (±15 %)
   ja taustatuulta (±15°) vaihdellaan satunnaisesti.
4. **Tulos** per CME: osumatodennäköisyys, saapumisen mediaani ja 5/10/90/95 %:n rajat, nopeus perillä.

Ensemble ajetaan vain, kun CME-lista muuttuu tai edellisestä ajosta on yli 12 h. Muulloin työ päättyy heti.

## Testaus paikallisesti

```bash
pip install -r <(grep -v '^huxt' requirements.txt) && pip install --no-deps huxt==5.0.3
HUXT_MEMBERS=50 python forecast.py --offline --now 2026-10-12T18:00:00Z --force
```

`--offline` käyttää kansion `testdata/` aineistoa: DONKI-simulaatiot 12.10.2026 ja Readingin
taustatuuli digitoituna heidän ennustekuvastaan. Tällä aineistolla tulosta voi verrata
Readingin julkaisemaan ennusteeseen (CME 618 km/s: osuma 65 %, saapuminen 14.10. 15:19 UTC).

## JSON-rakenne (lyhennetty)

```json
{
  "generated": "2026-10-12T18:00:00Z",
  "members": 500,
  "ambient_source": ["..."],
  "cmes": [
    { "id": "2026-10-11T23:12:00-CME-001", "v": 618, "lon": -28, "lat": -1, "width": 72,
      "status": "tulossa", "hit_pct": 65.0,
      "arrival_median": "2026-10-14T15:00Z", "arrival_p10": "...", "arrival_p90": "...",
      "speed_at_earth_median": 666 }
  ],
  "earth_speed_deterministic": [["2026-10-12T18:00Z", 640], ["...", 0]]
}
```

`status`: `tulossa`, `saapunut / ohi`, `ei osumaa`, `ei kohti Maata` tai `jo ohi (ei ajettu)`.

## Lähteet ja lisenssit

- HUXt: Barnard & Owens (2022), Frontiers in Physics, doi:10.3389/fphy.2022.1005621
- NASA DONKI (CCMC), NOAA SWPC: julkista dataa
- Tämä on kokeellinen harrasteennuste, ei virallinen avaruussääennuste.
