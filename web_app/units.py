UNITS = {
    "temperature": {
        "si": "K",
        "options": {
            "K": {"label": "K", "name": "kelvin", "to_si": lambda v: v, "from_si": lambda v: v, "decimals": 3},
            "C": {"label": "°C", "name": "grado Celsius", "to_si": lambda v: v + 273.15, "from_si": lambda v: v - 273.15, "decimals": 3},
        },
    },
    "pressure": {
        "si": "Pa",
        "options": {
            "Pa": {"label": "Pa", "name": "pascal", "to_si": lambda v: v, "from_si": lambda v: v, "decimals": 1},
            "kPa": {"label": "kPa", "name": "kilopascal", "to_si": lambda v: v * 1e3, "from_si": lambda v: v / 1e3, "decimals": 4},
            "bar": {"label": "bar", "name": "bar", "to_si": lambda v: v * 1e5, "from_si": lambda v: v / 1e5, "decimals": 5},
            "MPa": {"label": "MPa", "name": "megapascal", "to_si": lambda v: v * 1e6, "from_si": lambda v: v / 1e6, "decimals": 6},
        },
    },
    "density": {
        "si": "kg/m3",
        "options": {
            "kg/m3": {"label": "kg/m³", "name": "chilogrammo al metro cubo", "to_si": lambda v: v, "from_si": lambda v: v, "decimals": 4},
            "g/cm3": {"label": "g/cm³", "name": "grammo al centimetro cubo", "to_si": lambda v: v * 1e3, "from_si": lambda v: v / 1e3, "decimals": 7},
        },
    },
    "specific_energy": {
        "si": "J/kg",
        "options": {
            "J/kg": {"label": "J/kg", "name": "joule al chilogrammo", "to_si": lambda v: v, "from_si": lambda v: v, "decimals": 2},
            "kJ/kg": {"label": "kJ/kg", "name": "kilojoule al chilogrammo", "to_si": lambda v: v * 1e3, "from_si": lambda v: v / 1e3, "decimals": 5},
        },
    },
    "specific_entropy": {
        "si": "J/(kg*K)",
        "options": {
            "J/(kg*K)": {"label": "J/kg·K", "name": "joule al chilogrammo kelvin", "to_si": lambda v: v, "from_si": lambda v: v, "decimals": 3},
            "kJ/(kg*K)": {"label": "kJ/kg·K", "name": "kilojoule al chilogrammo kelvin", "to_si": lambda v: v * 1e3, "from_si": lambda v: v / 1e3, "decimals": 6},
        },
    },
    "dimensionless": {
        "si": "—",
        "options": {"—": {"label": "—", "name": "adimensionale", "to_si": lambda v: v, "from_si": lambda v: v, "decimals": 5}},
    },
    "speed": {
        "si": "m/s",
        "options": {"m/s": {"label": "m/s", "name": "metro al secondo", "to_si": lambda v: v, "from_si": lambda v: v, "decimals": 3}},
    },
    "inverse_temperature": {
        "si": "1/K",
        "options": {"1/K": {"label": "1/K", "name": "kelvin alla meno uno", "to_si": lambda v: v, "from_si": lambda v: v, "decimals": 8}},
    },
}

DEFAULTS = {
    "temperature": "K",
    "pressure": "bar",
    "density": "kg/m3",
    "specific_energy": "kJ/kg",
    "specific_entropy": "kJ/(kg*K)",
    "dimensionless": "—",
    "speed": "m/s",
    "inverse_temperature": "1/K",
}

# Units built from the chosen pressure and density: each kind gives the power of every base
# quantity, so the pressure menu also sets K/bar, 1/bar and bar/K and the density menu sets
# m³/kg or cm³/g. Temperature differences stay in kelvin whatever scale the input uses.
COMPOUND = {
    "inverse_pressure": {"pressure": -1},
    "temperature_per_pressure": {"kelvin": 1, "pressure": -1},
    "pressure_per_temperature": {"pressure": 1, "kelvin": -1},
    "pressure_per_density": {"pressure": 1, "density": -1},
    "pressure_per_density_squared": {"pressure": 1, "density": -2},
    "inverse_density": {"density": -1},
    "inverse_density_squared": {"density": -2},
}

_SYMBOLS = {
    "Pa": {"Pa": 1}, "kPa": {"kPa": 1}, "bar": {"bar": 1}, "MPa": {"MPa": 1},
    "kg/m3": {"kg": 1, "m": -3}, "g/cm3": {"g": 1, "cm": -3},
    "K": {"K": 1},
}

_SUPERSCRIPT = str.maketrans("-0123456789", "⁻⁰¹²³⁴⁵⁶⁷⁸⁹")


def spec(kind: str, unit: str) -> dict:
    group = UNITS[kind]
    try:
        chosen = group["options"].get(unit)
    except TypeError:
        chosen = None
    return chosen or group["options"][DEFAULTS[kind]]


def to_si(kind: str, unit: str, value):
    if value is None:
        return None
    return spec(kind, unit)["to_si"](value)


def from_si(kind: str, unit: str, value):
    if value is None:
        return None
    return spec(kind, unit)["from_si"](value)


def label(kind: str, unit: str) -> str:
    return spec(kind, unit)["label"]


def decimals(kind: str, unit: str) -> int:
    return spec(kind, unit)["decimals"]


def options(kind: str) -> list:
    return [{"label": o["label"], "value": key} for key, o in UNITS[kind]["options"].items()]


def _chosen(kind: str, display: dict) -> str:
    unit = (display or {}).get(kind)
    try:
        return unit if unit in UNITS[kind]["options"] else DEFAULTS[kind]
    except TypeError:
        return DEFAULTS[kind]


def _compound(kind: str, display: dict):
    powers, factor = {}, 1.0
    for base, power in COMPOUND[kind].items():
        if base == "kelvin":
            unit = "K"
        else:
            unit = _chosen(base, display)
            convert = spec(base, unit)
            factor *= (convert["from_si"](1.0) ** power if power > 0
                       else convert["to_si"](1.0) ** -power)
        for symbol, exponent in _SYMBOLS[unit].items():
            powers[symbol] = powers.get(symbol, 0) + exponent * power

    def written(symbol, power):
        return symbol if power == 1 else symbol + str(power).translate(_SUPERSCRIPT)

    above = [written(s, p) for s, p in powers.items() if p > 0]
    below = [written(s, -p) for s, p in powers.items() if p < 0]
    text = "·".join(above) or "1"
    return (f"{text}/{'·'.join(below)}" if below else text), factor


def display_label(kind: str, display: dict) -> str:
    if kind in COMPOUND:
        return _compound(kind, display)[0]
    return label(kind, _chosen(kind, display))


def display_value(kind: str, display: dict, value):
    if value is None:
        return None
    if kind in COMPOUND:
        return value * _compound(kind, display)[1]
    return from_si(kind, _chosen(kind, display), value)
