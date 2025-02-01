def test_import():
    import wickd


def test_orbital_space():
    import wickd

    wickd.reset_space()
    assert wickd.num_spaces() == 0


if __name__ == "__main__":
    test_import()
    test_orbital_space()
