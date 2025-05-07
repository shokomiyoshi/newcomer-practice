def hello_test1(input1: str, input2: str, switch1: bool, switch2: bool) -> None:
    """uvパッケージングのデモ用挨拶関数。

    指定されたスイッチの状態に応じて、入力された名前で挨拶メッセージを表示します。

    Args:
        input1: 最初の挨拶で使用する名前
        input2: 2番目の挨拶で使用する名前
        switch1: Trueの場合、input1を使った挨拶を表示
        switch2: Trueの場合、input2を使った挨拶を表示

    Returns:
        None

    Example:
        >>> hello_test1("太郎", "花子", True, False)
        Hello uv packaging!. I'm 太郎.
        >>> hello_test1("太郎", "花子", False, True)
        I'm 花子. Nice to meet you.
    """
    if switch1:
        print(f"Hello uv packaging!. I'm {input1}.")
    if switch2:
        print(f"I'm {input2}. Nice to meet you.")