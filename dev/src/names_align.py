from Bio.pairwise2 import align


class NamesAlign:
    @staticmethod
    def convert_to_chinese(a: list[str], b: list[str]):
        charcode = ord("一") - 1
        seqs = {seq: chr(charcode + ix) for ix, seq in enumerate(set(a + b))}
        return (
            "".join([seqs[seq] for seq in a]),
            "".join([seqs[seq] for seq in b]),
            seqs,
        )

    @staticmethod
    def chinese_to_list(a: str, seqs: dict[str, str]):
        seqs_rev = {v: k for k, v in seqs.items()}
        return [seqs_rev[seq] for seq in a]

    @classmethod
    def align(cls, new_names: list[str], base_names: list[str]):
        if not all((new_names, base_names)):
            return new_names or base_names
        a, b, seqs = cls.convert_to_chinese(new_names, base_names)
        alignments = align.globalxx(a, b)
        new_aligned, base_aligned = alignments[0][0], alignments[0][1]
        combined_aligned = [
            base_char if base_char != "-" else new_char
            for new_char, base_char in zip(new_aligned, base_aligned)
        ]
        combined = cls.chinese_to_list(combined_aligned, seqs)
        return combined
