from PIL import Image, ImageDraw, ImageFont
import argparse


def text_to_image(txt_file_path,
                  font_path='msyh.ttc',
                  font_size=50,
                  text_color='black',
                  bg_color='white'):
    # 读取文本文件内容
    with open(txt_file_path, 'r', encoding='utf-8') as f:
        text = f.read()
        text = text.replace('\\\\r\\\\n', '\\n')
        print(text)
    # 创建一个白色背景的空白图像
    img = Image.new('RGB', (1800, 1000), color=bg_color)

    # 在图像上创建一个Draw对象
    draw = ImageDraw.Draw(img)

    # 设置要绘制的文本和字体
    font = ImageFont.truetype(font_path, size=font_size)

    # 确定文本的位置，并使用指定的字体和颜色将其绘制到图像上
    # text_width, text_height = draw.textsize(text, font=font)
    textbox = draw.textbbox((0, 0), text, font=font)
    text_width = textbox[2] - textbox[0]
    text_height = textbox[3] - textbox[1]
    x = (img.width - text_width) // 2
    y = (img.height - text_height) // 2
    draw.text((x, y), text, fill=text_color, font=font)

    # 保存图像到文件并返回文件路径
    # 保存到figsave目录下，文件名为1_1.png, 1_2.png等
    txt_file_path = txt_file_path.replace('output_part_', '1_')
    output_file_path = 'figsave/' + txt_file_path.replace('.txt', '.jpg')
    img.save(output_file_path)
    print('success')
    return output_file_path


# 将output_part_1.txt, output_part_2.txt等文件转换为图片1_1.jpg, 1_2.jpg等
for i in range(4):
    text_to_image(f'output_part_{i+1}.txt',
                  font_path='usr/share/fonts/opentype/noto/NotoSansCJK-Regular.ttc',
                  font_size=50,
                  text_color='black',
                  bg_color='white')
