import pytest
from chemicalgof import decode, encode, split_fragsmiles

def _extract_stereo_labels(sampled:list[str]):
    stereo_elements :list[str] = []
    for element in sampled:
        if (
            element.startswith('<') and ('R' in element.upper() or 'S' in element.upper())
        ) or (
            '|' in element
        ):
            stereo_elements.append(element)

    return tuple(stereo_elements)

@pytest.mark.parametrize(
    "sampled",
    [

        [
            "C",
            "<12R>",
            "O=C1NCCCc2cccc(c2)CCCOCc2cccc1c2",
            "<20R>",
            "<16>",
            "(",
            "C",
            ")",
            "<11R>",
            "(",
            "<6>",
            "c1ccc2c(c1)CCC21CCNCC1",
            ")",
            "<22S>",
            "(",
            "O",
            ")",
            "O",
            "C",
        ],
    ],
    ids=[
        "required_unstricted_chirality",
    ],
)

def test_decoding_sampled(sampled):

    with pytest.raises(Exception) as exception_info:
        decoded_stricted = decode(sampled, strict_chirality=True)

    assert "Chirality Error" in str(exception_info.value)

    decoded_unstricted = decode(sampled, strict_chirality=False)
    reencoded = encode(decoded_unstricted)
    reencoded_splitted = split_fragsmiles(reencoded)

    requested_stereo_elements = _extract_stereo_labels(sampled)
    actual_stereo_elements = _extract_stereo_labels(reencoded_splitted)

    assert len(requested_stereo_elements) != len(actual_stereo_elements) or requested_stereo_elements != actual_stereo_elements
