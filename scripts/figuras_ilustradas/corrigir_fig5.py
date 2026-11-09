import argparse
from pathlib import Path

from PIL import Image, ImageDraw, ImageFont

PAINEIS = [(0, 632), (640, 1278), (1288, 1920)]
PONTOS = [(-0.22, -0.22), (-0.19, 0.0), (-0.91, -0.63)]
Y0, Y1 = 760, 965
CINZA, PONTO, LAR = (60, 60, 60), (110, 115, 130), (230, 120, 60)


def fonte(tam=34):
    for nome in ('arial.ttf', 'DejaVuSans.ttf'):
        try:
            return ImageFont.truetype(nome, tam)
        except OSError:
            pass
    return ImageFont.load_default()


def corrigir(entrada, saida):
    im = Image.open(entrada).convert('RGB')
    d, f = ImageDraw.Draw(im), fonte()
    for (x0, x1), (sem, com) in zip(PAINEIS, PONTOS):
        d.rectangle([x0 + 4, Y0, x1 - 4, Y1], fill=im.getpixel((x0 + 20, 735)))
        xz, escala = x0 + int((x1 - x0) * 0.72), (x1 - x0) * 0.30
        d.line([xz, Y0 + 30, xz, Y1 - 15], fill=CINZA, width=3)
        d.text((xz - 6, Y1 - 12), '0', fill=CINZA, font=f, anchor='rt')
        for y, rot, v in [(Y0 + 75, 'sem ajuste', sem), (Y0 + 140, 'com ajuste', com)]:
            d.text((x0 + 30, y), rot, fill=(30, 30, 30), font=f, anchor='lm')
            cx = xz + int(v * escala)
            d.ellipse([cx - 14, y - 14, cx + 14, y + 14], fill=PONTO)
        if com != 0 and abs(com) < abs(sem) - 0.1:
            cx = xz + int(com * escala) + 34
            d.polygon([(cx, Y0 + 155), (cx + 30, Y0 + 155), (cx + 15, Y0 + 125)], fill=LAR)
    im.save(saida)


if __name__ == '__main__':
    ap = argparse.ArgumentParser(description='Redesenha a faixa inferior da ilustração das duas metades com os g do artigo.')
    ap.add_argument('entrada', type=Path)
    ap.add_argument('saida', type=Path)
    a = ap.parse_args()
    corrigir(a.entrada, a.saida)
