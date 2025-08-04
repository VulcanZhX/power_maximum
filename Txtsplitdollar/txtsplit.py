def split_text_file_by_dollar(input_file, output_prefix):
    with open(input_file, 'r', encoding='utf-8') as f:
        content = f.read()

    # 按 $ 分割文本
    parts = content.split('$')

    # 去除每部分前后空白，并过滤空字符串
    parts = [part.strip() for part in parts if part.strip()]

    # 将每部分写入单独文件
    for i, part in enumerate(parts, 1):
        output_file = f"{output_prefix}_{i}.txt"
        with open(output_file, 'w', encoding='utf-8') as f_out:
            f_out.write(part)
        print(f"Saved part {i} to {output_file}")


# 示例用法
if __name__ == "__main__":
    input_path = "input.txt"  # 你的输入文件名
    output_prefix = "output_part"  # 输出文件名前缀
    split_text_file_by_dollar(input_path, output_prefix)
