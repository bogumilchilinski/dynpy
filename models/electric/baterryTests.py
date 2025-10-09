from __future__ import annotations

import argparse
import re
import sys
import time
from datetime import datetime
from pathlib import Path
from typing import Optional, Tuple


try:
    import serial
    import serial.tools.list_ports
except Exception:
    serial = None  


try:
    import pandas as pd
except Exception:
    pd = None

try:
    import matplotlib.pyplot as plt
except Exception:
    plt = None


# === code generators ===
def get_pico_script() -> str:

    code = """

      from machine import ADC, Pin
import time
adc = ADC(Pin(26))  
VREF = 3.26
MAX_ADC = 65535
R1 = 9780
R2 = 9790
K = 1.0675
def adc_avg(n=16):
    return sum(adc.read_u16() for _ in range(n)) 
print("Kalibracja offsetu (5 sek)...")
offset_sum = 0
samples = 0
start = time.ticks_ms()
while time.ticks_diff(time.ticks_ms(), start) < 5000:
    offset_sum += adc_avg()
    samples += 1
    time.sleep(0.01)
offset = offset_sum / samples
print("Zakonczono Offset = {:.0f}".format(offset))
def odczytaj_napiecie():
    raw = adc_avg()
    raw_corr = max(0, raw - offset)
    v_adc = (raw_corr / MAX_ADC) * VREF
    v_bat = v_adc * (R1 + R2) / R2
    v_bat *= K  
    return v_bat

try:
    while True:
        napiecie = odczytaj_napiecie()
        print("Napiecie: {:.3f} V".format(napiecie))
        time.sleep(1)

except KeyboardInterrupt:
    print("Zatrzymano")

    """
    return code


def get_host_script() -> str:
    code = """import serial
import serial.tools.list_ports
import pandas as pd
import matplotlib.pyplot as plt
from datetime import datetime
import time


def find_pico_port():
    ports = serial.tools.list_ports.comports()
    for port in ports:
        if "ttyACM" in port.device or "usbmodem" in port.device:
            return port.device
    raise RuntimeError("Nie znaleziono Pico. Upewnij się, że jest podłączone przez USB.")

def calculate_soc(voltage):
    if voltage >= 4.20:
        return 1.00
    elif voltage >= 4.18:
        return 0.98
    elif voltage >= 4.16:
        return 0.96
    elif voltage >= 4.14:
        return 0.94
    elif voltage >= 4.12:
        return 0.92
    elif voltage >= 4.10:
        return 0.90
    elif voltage >= 4.08:
        return 0.88
    elif voltage >= 4.06:
        return 0.86
    elif voltage >= 4.04:
        return 0.84
    elif voltage >= 4.02:
        return 0.82
    elif voltage >= 4.00:
        return 0.80
    elif voltage >= 3.98:
        return 0.78
    elif voltage >= 3.96:
        return 0.76
    elif voltage >= 3.94:
        return 0.74
    elif voltage >= 3.92:
        return 0.72
    elif voltage >= 3.90:
        return 0.70
    elif voltage >= 3.88:
        return 0.68
    elif voltage >= 3.86:
        return 0.66
    elif voltage >= 3.84:
        return 0.64
    elif voltage >= 3.82:
        return 0.62
    elif voltage >= 3.80:
        return 0.60
    elif voltage >= 3.78:
        return 0.58
    elif voltage >= 3.76:
        return 0.56
    elif voltage >= 3.74:
        return 0.54
    elif voltage >= 3.72:
        return 0.52
    elif voltage >= 3.70:
        return 0.50
    elif voltage >= 3.68:
        return 0.47
    elif voltage >= 3.66:
        return 0.44
    elif voltage >= 3.64:
        return 0.41
    elif voltage >= 3.62:
        return 0.38
    elif voltage >= 3.60:
        return 0.35
    elif voltage >= 3.58:
        return 0.32
    elif voltage >= 3.56:
        return 0.29
    elif voltage >= 3.54:
        return 0.26
    elif voltage >= 3.52:
        return 0.23
    elif voltage >= 3.50:
        return 0.20
    elif voltage >= 3.48:
        return 0.17
    elif voltage >= 3.45:
        return 0.14
    elif voltage >= 3.42:
        return 0.11
    elif voltage >= 3.40:
        return 0.08
    elif voltage >= 3.38:
        return 0.06
    elif voltage >= 3.36:
        return 0.04
    elif voltage >= 3.34:
        return 0.02
    else:
        return 0.0





try:
    pico_port = find_pico_port()
    print(f"[INFO] Wykryto Pico na porcie: {pico_port}")
except RuntimeError as e:
    print(f"[ERROR] {e}")
    exit(1)


BAUDRATE = 115200
timestamp = datetime.now().strftime('%Y-%m-%d_%H-%M-%S')
csv_file = f'pomiary_{timestamp}.csv'
excel_file = f'pomiary_{timestamp}.xlsx'
plot_soc_time = f'soc_vs_time_{timestamp}.png'
plot_soc_voltage = f'soc_vs_voltage_{timestamp}.png'
data = []

# === Otwórz połączenie szeregowe ===
try:
    ser = serial.Serial(pico_port, BAUDRATE, timeout=1)
    time.sleep(2)
    print("[INFO] Rozpoczynam odczyt danych z Pico (Ctrl+C aby zakończyć)...")

    while True:
        line = ser.readline().decode("utf-8").strip()
        if line.startswith("Napięcie:"):
            try:
                voltage = float(line.split()[1])
                soc = calculate_soc(voltage)
                now = datetime.now()
                print(f"{now.strftime('%H:%M:%S')} - {voltage:.3f} V, SOC: {soc*100:.1f}%")
                data.append((now, voltage, soc))
            except:
                continue

except KeyboardInterrupt:
    print("\n[INFO] Przerwano. Zapisuję dane...")

    df = pd.DataFrame(data, columns=["czas", "napięcie", "SOC"])
    df.to_csv(csv_file, index=False)
    df.to_excel(excel_file, index=False)


    plt.figure()
    plt.plot(df["czas"], df["SOC"])
    plt.xlabel("Czas")
    plt.ylabel("SOC (0–1)")
    plt.title("SOC od czasu")
    plt.xticks(rotation=45)
    plt.grid()
    plt.tight_layout()
    plt.savefig(plot_soc_time)


    plt.figure()
    plt.plot(df["napięcie"], df["SOC"])
    plt.xlabel("Napięcie [V]")
    plt.ylabel("SOC (0–1)")
    plt.title("SOC od napięcia")
    plt.grid()
    plt.tight_layout()
    plt.savefig(plot_soc_voltage)

    print(f"[INFO] Dane zapisane do: {csv_file}, {excel_file}")
    print(f"[INFO] Wykresy: {plot_soc_time}, {plot_soc_voltage}")
"""
    return code



def _require(dep, name: str):
    if dep is None:
        raise RuntimeError(
            f"Ta funkcja wymaga pakietu '{name}'. Zainstaluj: pip install {name}"
        )


def find_pico_port() -> str:

    _require(serial, "pyserial")
    ports = serial.tools.list_ports.comports()
    for p in ports:
        dev = (p.device or "").lower()
        desc = (p.description or "").lower()
        if "ttyacm" in dev or "usbmodem" in dev or "raspberry" in desc or "pico" in desc:
            return p.device
    # Fallback: pierwszy port wyglądający na USB/COM
    for p in ports:
        dev = (p.device or "").lower()
        if dev.startswith("com") or "usb" in dev:
            return p.device
    raise RuntimeError("Nie znaleziono Pico. Upewnij się, że jest podłączone przez USB.")


def _parse_voltage(line: str) -> Optional[float]:

    # Szybka ścieżka: spróbuj wyłuskać liczbę po słowie
    m = re.search(r"([-+]?\d+(?:\.\d+)?)", line.replace(",", "."))
    if not m:
        return None
    try:
        return float(m.group(1))
    except ValueError:
        return None


def calculate_soc_default(voltage: float) -> float:

    thresholds = [
        (4.20, 1.00), (4.18, 0.98), (4.16, 0.96), (4.14, 0.94), (4.12, 0.92),
        (4.10, 0.90), (4.08, 0.88), (4.06, 0.86), (4.04, 0.84), (4.02, 0.82),
        (4.00, 0.80), (3.98, 0.78), (3.96, 0.76), (3.94, 0.74), (3.92, 0.72),
        (3.90, 0.70), (3.88, 0.68), (3.86, 0.66), (3.84, 0.64), (3.82, 0.62),
        (3.80, 0.60), (3.78, 0.58), (3.76, 0.56), (3.74, 0.54), (3.72, 0.52),
        (3.70, 0.50), (3.68, 0.47), (3.66, 0.44), (3.64, 0.41), (3.62, 0.38),
        (3.60, 0.35), (3.58, 0.32), (3.56, 0.29), (3.54, 0.26), (3.52, 0.23),
        (3.50, 0.20), (3.48, 0.17), (3.45, 0.14), (3.42, 0.11), (3.40, 0.08),
        (3.38, 0.06), (3.36, 0.04), (3.34, 0.02),
    ]
    for thr, soc in thresholds:
        if voltage >= thr:
            return soc
    return 0.0


def collect_data(
    baudrate: int = 115200,
    out_dir: Path | str = ".",
    soc_fn = calculate_soc_default,
    save_excel: bool = True,
    make_plots: bool = True,
) -> Tuple[Path, Optional[Path], Optional[Path], Optional[Path]]:

    _require(serial, "pyserial")
    if save_excel:
        _require(pd, "pandas")
    if make_plots:
        _require(plt, "matplotlib")

    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    pico_port = find_pico_port()
    print(f"[INFO] Wykryto Pico na porcie: {pico_port}")

    timestamp = datetime.now().strftime("%Y-%m-%d_%H-%M-%S")
    csv_path = out_dir / f"pomiary_{timestamp}.csv"
    xlsx_path = out_dir / f"pomiary_{timestamp}.xlsx"
    plot_soc_time = out_dir / f"soc_vs_time_{timestamp}.png"
    plot_soc_voltage = out_dir / f"soc_vs_voltage_{timestamp}.png"

    ser = serial.Serial(pico_port, baudrate, timeout=1)
    time.sleep(2)
    print("[INFO] Rozpoczynam odczyt danych (Ctrl+C aby zakończyć)...")

    rows = []
    try:
        while True:
            raw = ser.readline()
            if not raw:
                continue
            try:
                line = raw.decode("utf-8", errors="ignore").strip()
            except Exception:
                continue

            v = _parse_voltage(line)
            if v is None:
                continue

            soc = float(soc_fn(v))
            now = datetime.now()
            print(f"{now.strftime('%H:%M:%S')} - {v:.3f} V, SOC: {soc*100:.1f}%")
            rows.append((now.isoformat(), v, soc))

    except KeyboardInterrupt:
        print("\n[INFO] Przerwano. Zapisuję dane...")

    finally:
        try:
            ser.close()
        except Exception:
            pass


    with csv_path.open("w", encoding="utf-8") as f:
        f.write("czas,napięcie,SOC\n")
        for t, v, s in rows:
            f.write(f"{t},{v:.6f},{s:.6f}\n")


    if save_excel and pd is not None:
        df = pd.DataFrame(rows, columns=["czas", "napięcie", "SOC"])
        df.to_excel(xlsx_path, index=False)
    else:
        xlsx_path = None


    if make_plots and plt is not None:
        # SOC vs czas
        ts = [datetime.fromisoformat(t) for (t, _, _) in rows]
        socs = [s for (_, _, s) in rows]
        plt.figure()
        plt.plot(ts, socs)
        plt.xlabel("Czas")
        plt.ylabel("SOC (0–1)")
        plt.title("SOC od czasu")
        plt.xticks(rotation=45)
        plt.grid()
        plt.tight_layout()
        plt.savefig(plot_soc_time)

        # SOC vs napięcie
        volts = [v for (_, v, _) in rows]
        plt.figure()
        plt.plot(volts, socs)
        plt.xlabel("Napięcie [V]")
        plt.ylabel("SOC (0–1)")
        plt.title("SOC od napięcia")
        plt.grid()
        plt.tight_layout()
        plt.savefig(plot_soc_voltage)
    else:
        plot_soc_time = None
        plot_soc_voltage = None

    print(f"[INFO] CSV:   {csv_path}")
    if xlsx_path:
        print(f"[INFO] XLSX:  {xlsx_path}")
    if plot_soc_time:
        print(f"[INFO] Wykresy: {plot_soc_time}, {plot_soc_voltage}")

    return csv_path, xlsx_path, plot_soc_time, plot_soc_voltage



def _build_cli() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(description="Battery test helper")
    sub = p.add_subparsers(dest="cmd", required=True)

    sub.add_parser("print-pico", help="Wypisz kod MicroPython na Pico")
    sub.add_parser("print-host", help="Wypisz kod hosta (RPi/PC)")

    c = sub.add_parser("collect", help="Uruchom kolektor danych (modułowy)")
    c.add_argument("--baudrate", type=int, default=115200)
    c.add_argument("--out-dir", type=str, default=".")
    c.add_argument("--no-excel", action="store_true", help="Nie zapisuj XLSX")
    c.add_argument("--no-plots", action="store_true", help="Nie generuj wykresów")

    return p


def main(argv=None) -> int:
    args = _build_cli().parse_args(argv)

    if args.cmd == "print-pico":
        print(get_pico_script())
        return 0

    if args.cmd == "print-host":
        print(get_host_script())
        return 0

    if args.cmd == "collect":
        try:
            collect_data(
                baudrate=getattr(args, "baudrate", 115200),
                out_dir=getattr(args, "out_dir", "."),
                save_excel=not getattr(args, "no_excel", False),
                make_plots=not getattr(args, "no_plots", False),
            )
            return 0
        except Exception as e:
            print(f"[ERROR] {e}")
            return 2

    return 0


if __name__ == "__main__":
    raise SystemExit(main())



from dynpy.utilities.report import Picture, NoEscape

def get_battery_test_tikz() -> str:

    return r"""
\begin{circuitikz}[american]

  \node[draw, rounded corners=6pt, minimum width=3.2cm, minimum height=1.2cm] (rpi4) at (0,0) {RPi 4};
  \node[draw, rounded corners=6pt, minimum width=3.8cm, minimum height=1.2cm, right=3.2cm of rpi4] (pico) {RPI PICO};
  \node[draw, rounded corners=6pt, minimum width=3.2cm, minimum height=1.2cm, right=5.5cm of pico] (motor) {DC Motor};

  \draw (rpi4.east) -- ($(rpi4.east)!0.5!(pico.west)$) node[midway, above,yshift=2pt]{UART} -- (pico.west);

  \coordinate (gp26) at ($(pico.east) + (1.0,0)$);
  \draw (pico.east) -- (gp26);
  \node[above] at (gp26) {\small GP26};

  \draw (gp26) to[R,l=$R_{1}$] ++(0,-2.2) node[ground]{};
  \draw (gp26) to[R,l=$R_{2}$] ++(2.2,0) coordinate (afterR2);
  \draw (afterR2) -- (motor.west);

  \coordinate (toSwitch) at ($(motor.east) + (1.0,0)$);
  \draw (motor.east) -- (toSwitch)
        to[spst,l^=$S$] ++(2.0,0) coordinate (rightBus);

  \draw (rightBus) -- ++(0,-2.5) coordinate (batTop);
  \draw (batTop) to[battery1,l_=18650] ++(0,-2.2) node[ground]{};

  \draw (motor.east) node[circ]{};
  \draw (rightBus) node[circ]{};
\end{circuitikz}
"""


def get_battery_test_picture() -> Picture:

    return Picture(
        fr'dynpy\models\images\baterryTestElectricScheme.png',
        caption='Schemat elektryczny stanowiska testowego (RPi4 – Pico – DC Motor – 18650)',
        width=NoEscape(r'0.65\textwidth')
    )




# display(get_battery_test_picture())


# display(TikZPicture(get_battery_test_tikz()))