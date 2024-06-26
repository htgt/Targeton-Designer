from dataclasses import dataclass


@dataclass
class Interval:
    start: int
    end: int


@dataclass
class DesignedPrimer:
    name: str
    penalty: float
    pair_id: str
    sequence: str
    coords: Interval
    primer_start: int
    primer_end: int
    strand: str
    tm: float
    gc_percent: float
    self_any_th: float
    self_end_th: float
    hairpin_th: float
    end_stability: float

    def __eq__(self, other):
        if isinstance(other, DesignedPrimer):
            # Exclude the 'name' and 'pair_id' attributes from the comparison
            return (
                    self.penalty == other.penalty and
                    self.sequence == other.sequence and
                    self.coords == other.coords and
                    self.primer_start == other.primer_start and
                    self.primer_end == other.primer_end and
                    self.strand == other.strand and
                    self.tm == other.tm and
                    self.gc_percent == other.gc_percent and
                    self.self_any_th == other.self_any_th and
                    self.self_end_th == other.self_end_th and
                    self.hairpin_th == other.hairpin_th and
                    self.end_stability == other.end_stability
            )
        return False

    def __hash__(self):
        # Exclude the 'name' and 'pair_id' attributes from the comparison
        return hash((
            self.penalty,
            self.sequence,
            (self.coords.start, self.coords.end),
            self.primer_start,
            self.primer_end,
            self.strand,
            self.tm,
            self.gc_percent,
            self.self_any_th,
            self.self_end_th,
            self.hairpin_th,
            self.end_stability
        ))


def map_to_designed_primer(primer: dict):
    return DesignedPrimer(
        name=primer["primer"],
        penalty=primer["penalty"],
        pair_id=primer["pair_id"],
        sequence=primer["sequence"],
        coords=Interval(primer["coords"][0], primer["coords"][1]),
        primer_start=primer["primer_start"],
        primer_end=primer["primer_end"],
        strand=primer["strand"],
        tm=primer["tm"],
        gc_percent=primer["gc_percent"],
        self_any_th=primer["self_any_th"],
        self_end_th=primer["self_end_th"],
        hairpin_th=primer["hairpin_th"],
        end_stability=primer["end_stability"]
    )
