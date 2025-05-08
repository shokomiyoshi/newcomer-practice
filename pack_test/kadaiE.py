from Bio import SeqIO
import gzip
import argparse

# 残基数がcutoff以下の割合を計算して表示する関数
def calculate_ratio(fastafile: str, cutoff: int) -> None:

    count = 0
    allrecords = 0

    # ファイルを開いてFASTA形式で読み込む
    with gzip.open(fastafile, "rt") as handle:
        for record in SeqIO.parse(handle, 'fasta'):
            # ヘッダー行から情報を抽出
            header = record.description

            # mol:protein が含まれている場合
            if "mol:protein" in header:
                # "length:" を基に長さの値を抽出
                length_str = header.split("length:")[1].split()[0]  # length: の後の値を抽出
                length = int(length_str)  # 数値に変換

                # cutoff以下でカウント
                if length <= cutoff:
                    count += 1

                allrecords += 1  # トータルのシーケンス数をカウント

    # 比率を計算
    if allrecords > 0:
        ratio_less_than_cutoff = count / allrecords * 100
        print(f"{cutoff}残基以下の割合: {ratio_less_than_cutoff:.2f}%")
    else:
        print("mol:protein を含むレコードがありませんでした。")


# CLI引数からファイルとカットオフ値を受け取る関数
def kadaie_1():
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input", required=True, type=str, help="pdb_seqres.txt.gz ファイルのパス")
    parser.add_argument("-l", "--length", required=True, type=int, help="残基数のカットオフ（例：100）")
    args = parser.parse_args()
    calculate_ratio(args.input, args.length)

if __name__ == "__main__":
    kadaie_1()