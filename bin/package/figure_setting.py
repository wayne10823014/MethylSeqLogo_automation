import math
import matplotlib
matplotlib.use('Agg')  # 不显示图形界面
import matplotlib.pyplot as plt
import matplotlib.patheffects as patheffects
from matplotlib.transforms import Affine2D
from matplotlib.gridspec import GridSpec
from matplotlib.font_manager import FontProperties
import matplotlib.patches as patches
from matplotlib import transforms
import seaborn as sns
import numpy as np

# 定义颜色方案
COLOR_SCHEME = {
    'G': 'orange',
    'A': 'green',
    'C': 'deepskyblue',
    'T': 'red'
}

BAR_COLOR_SCHEME = {
    'CG': 'black',
    'CHG': 'gray',
    'CHH': 'lightgray',
}

def set_fig(entropys, logotype, mode, plotlen):
    """
    设置图形大小。

    参数：
    - entropys (DataFrame): 信息熵数据。
    - logotype (str): 图形类型。
    - mode (str): 模式。
    - plotlen (int): 绘图长度。

    返回：
    - fig (Figure): Matplotlib 图形对象。
    """
    if logotype == 'riverlake':
        Height = 3.0
        figureheight = 3.0
    else:
        if mode == 'Methyl':
            Height = 3.0
            figureheight = Height + 4.0
        else:
            Height = max(math.ceil(max(entropys['Base'])), 3.0)
            figureheight = Height
    fig = plt.figure(figsize=(plotlen + 1, figureheight))
    return fig

class Scale(patheffects.RendererBase):
    """
    自定义缩放效果，用于调整文本的大小。
    """
    def __init__(self, sx, sy=None):
        self._sx = sx
        self._sy = sy

    def draw_path(self, renderer, gc, tpath, affine, rgbFace):
        affine = Affine2D().scale(self._sx, self._sy) + affine
        renderer.draw_path(gc, tpath, affine, rgbFace)

class SeqLogoPlot:
    """
    序列标志图绘制类。
    """
    def __init__(
        self, fig, celltype, four_base_heights, entropys, Cmethyls, Gmethyls, bgpps,
        dientropys, bg_dientropys_max, bg_dientropys_min, J_bCG, J_bCHG, J_bCHH,
        Freqs_, mode, plotlen, threshold, TF
    ):
        self.fig = fig
        self.celltype = celltype
        self.base_heights = four_base_heights
        self.entropys = entropys
        self.Cmethyls = Cmethyls
        self.Gmethyls = Gmethyls
        self.bgpps = bgpps
        self.dientropys = dientropys
        self.bg_dientropys_max = bg_dientropys_max
        self.bg_dientropys_min = bg_dientropys_min
        self.J_bCG = J_bCG
        self.J_bCHG = J_bCHG
        self.J_bCHH = J_bCHH
        self.Freqs_ = Freqs_
        self.mode = mode
        self.plotlen = plotlen
        self.threshold = threshold
        self.TF = TF
        print("Threshold:", threshold)

    def plot_logo(self):
        """
        绘制序列标志图。
        """
        # 清除图形内容
        self.fig.clf()
        self.fig.tight_layout()

        if self.mode == 'Methyl':
            Height = 3.0
            figureheight = Height + 4.0
        else:
            Height = max(math.ceil(max(self.entropys['Base'])), 3.0)
            figureheight = Height

        if self.mode == 'Methyl':
            gs = GridSpec(
                3,  # 行数
                2,  # 列数
                width_ratios=[self.plotlen + 1, 1],
                height_ratios=[3, 6, 3],
                wspace=0.0,
                hspace=0.5,
            )

            ax1 = plt.subplot(gs[2])
            ax2 = plt.subplot(gs[3], sharey=ax1)
            ax3 = plt.subplot(gs[0], sharex=ax1)
            ax4 = plt.subplot(gs[1], sharex=ax2)
            ax5 = plt.subplot(gs[4], sharex=ax1)
            ax6 = plt.subplot(gs[5], sharey=ax5)

            x_axis = list(range(self.plotlen))
            ax1.set_xticks(x_axis)
            ax1.set_xticklabels([str(i + 1) for i in x_axis])

            y_axis = range(int(Height) + 1)

            font_T = FontProperties(size=68, weight='bold', family='monospace')
            font_C = FontProperties(size=66, weight='bold', family='monospace')

            bbox_props = dict(boxstyle="square, pad=0.0", fill=False, lw=0.0, alpha=0.5)

            base_heights = zip(self.base_heights, self.Cmethyls, self.Gmethyls)

            for i, (bases, Cmethyl, Gmethyl) in enumerate(base_heights):
                Bottom = 0
                width = 1

                for (base, score) in bases:
                    ax1.bar(
                        x_axis[i],
                        score,
                        bottom=Bottom,
                        width=width,
                        color=COLOR_SCHEME[base],
                        alpha=0.0
                    )

                    if base in ['C', 'G']:
                        yshift = 1.0
                        txt1 = ax1.text(
                            x_axis[i],
                            Bottom,
                            base,
                            fontproperties=font_C,
                            transform=transforms.offset_copy(
                                ax1.transData,
                                fig=self.fig,
                                x=-6,
                                y=yshift * score,
                                units='points'
                            ),
                            ha='center',
                            va='baseline',
                            color=COLOR_SCHEME[base],
                            bbox=bbox_props,
                            zorder=1,
                        )
                        txt2 = ax1.text(
                            x_axis[i],
                            Bottom,
                            base,
                            fontproperties=font_C,
                            transform=transforms.offset_copy(
                                ax1.transData,
                                fig=self.fig,
                                x=-6,
                                y=yshift * score,
                                units='points'
                            ),
                            clip_on=True,
                            ha='center',
                            va='baseline',
                            color='black',
                            bbox=bbox_props,
                            zorder=2,
                        )
                        for txt in [txt1, txt2]:
                            txt.set_path_effects([Scale(1.2, score)])

                        # 计算甲基化水平
                        if base == 'C':
                            bgMlevel = score * (
                                self.Freqs_.iloc[i]['CpG_p'] * self.J_bCG +
                                self.Freqs_.iloc[i]['CHG_p'] * self.J_bCHG +
                                self.Freqs_.iloc[i]['CHH_p'] * self.J_bCHH
                            )
                            foreMlevel = score * (
                                self.Freqs_.iloc[i]['CpG_p'] * Cmethyl[0] +
                                self.Freqs_.iloc[i]['CHG_p'] * Cmethyl[1] +
                                self.Freqs_.iloc[i]['CHH_p'] * Cmethyl[2]
                            )
                        elif base == 'G':
                            bgMlevel = score * (
                                self.Freqs_.iloc[i]['CpG_m'] * self.J_bCG +
                                self.Freqs_.iloc[i]['CHG_m'] * self.J_bCHG +
                                self.Freqs_.iloc[i]['CHH_m'] * self.J_bCHH
                            )
                            foreMlevel = score * (
                                self.Freqs_.iloc[i]['CpG_m'] * Gmethyl[0] +
                                self.Freqs_.iloc[i]['CHG_m'] * Gmethyl[1] +
                                self.Freqs_.iloc[i]['CHH_m'] * Gmethyl[2]
                            )

                        # 绘制甲基化阴影
                        bot = Bottom
                        left = x_axis[i] - 0.5

                        p1 = patches.Rectangle(
                            (left, bot),
                            1.0,
                            foreMlevel,
                            clip_on=True,
                            lw=0,
                            fill=False
                        )

                        ax1.add_patch(p1)
                        txt2.set_clip_path(p1)

                        if score >= self.threshold:
                            ax1.plot(
                                (left, left + 1),
                                (bot + bgMlevel, bot + bgMlevel),
                                linewidth=1.0,
                                linestyle='--',
                                color='lightgray',
                                zorder=3,
                            )
                    else:
                        txt1 = ax1.text(
                            x_axis[i],
                            Bottom,
                            base,
                            fontproperties=font_T,
                            ha='center',
                            va='baseline',
                            color=COLOR_SCHEME[base],
                            bbox=bbox_props
                        )
                        txt1.set_path_effects([Scale(1.0, score)])

                    Bottom += score

            # 添加图例到 ax2
            fontKey = FontProperties(size=19, weight='bold', family='monospace')
            ax2top = int(Height)
            methylkeyori = ax2top - 1.0
            ax2.set_xlim(0.0, 1.2)
            ax2.set_ylim(0.0, int(Height))

            # 绘制 CG 图例
            CGC = ax2.text(0.2, methylkeyori + 0.65, 'C', color='white', ha='center', fontproperties=fontKey,
                           transform=transforms.offset_copy(
                               ax2.transData,
                               fig=self.fig,
                               x=0.0,
                               y=1.5,
                               units='points'
                           ))
            CGG = ax2.text(0.5, methylkeyori + 0.65, 'G', color='white', ha='center', fontproperties=fontKey,
                           transform=transforms.offset_copy(
                               ax2.transData,
                               fig=self.fig,
                               x=0.0,
                               y=1.5,
                               units='points'
                           ))
            cgc = ax2.text(0.2, methylkeyori + 0.65, 'C', color='black', ha='center', fontproperties=fontKey,
                           clip_on=True,
                           transform=transforms.offset_copy(
                               ax2.transData,
                               fig=self.fig,
                               x=0.0,
                               y=1.5,
                               units='points'
                           ))
            cgg = ax2.text(0.5, methylkeyori + 0.65, 'G', color='black', ha='center', fontproperties=fontKey,
                           clip_on=True,
                           transform=transforms.offset_copy(
                               ax2.transData,
                               fig=self.fig,
                               x=0.0,
                               y=1.5,
                               units='points'
                           ))

            for txt in [CGC, CGG, cgc, cgg]:
                txt.set_path_effects([patheffects.Stroke(linewidth=2, foreground='black'), patheffects.Normal()])

            p1 = patches.Rectangle(
                (0.0, methylkeyori + 0.65),
                1.0,
                self.J_bCG * 0.33,
                clip_on=True,
                lw=0,
                fill=False
            )
            ax2.add_patch(p1)
            cgc.set_clip_path(p1)
            cgg.set_clip_path(p1)

            # 绘制 CHG 图例
            CHGC = ax2.text(0.2, methylkeyori + 0.325, 'C', color='white', ha='center', fontproperties=fontKey,
                            transform=transforms.offset_copy(
                                ax2.transData,
                                fig=self.fig,
                                x=0.0,
                                y=1.5,
                                units='points'
                            ))
            CHGH = ax2.text(0.5, methylkeyori + 0.325, 'H', color='white', ha='center', fontproperties=fontKey,
                            transform=transforms.offset_copy(
                                ax2.transData,
                                fig=self.fig,
                                x=0.0,
                                y=1.5,
                                units='points'
                            ))
            CHGG = ax2.text(0.8, methylkeyori + 0.325, 'G', color='white', ha='center', fontproperties=fontKey,
                            transform=transforms.offset_copy(
                                ax2.transData,
                                fig=self.fig,
                                x=0.0,
                                y=1.5,
                                units='points'
                            ))

            chgc = ax2.text(0.2, methylkeyori + 0.325, 'C', color='black', ha='center', fontproperties=fontKey,
                            clip_on=True,
                            transform=transforms.offset_copy(
                                ax2.transData,
                                fig=self.fig,
                                x=0.0,
                                y=1.5,
                                units='points'
                            ))
            chgh = ax2.text(0.5, methylkeyori + 0.325, 'H', color='black', ha='center', fontproperties=fontKey,
                            clip_on=True,
                            transform=transforms.offset_copy(
                                ax2.transData,
                                fig=self.fig,
                                x=0.0,
                                y=1.5,
                                units='points'
                            ))
            chgg = ax2.text(0.8, methylkeyori + 0.325, 'G', color='black', ha='center', fontproperties=fontKey,
                            clip_on=True,
                            transform=transforms.offset_copy(
                                ax2.transData,
                                fig=self.fig,
                                x=0.0,
                                y=1.5,
                                units='points'
                            ))

            for txt in [CHGC, CHGH, CHGG, chgc, chgh, chgg]:
                txt.set_path_effects([patheffects.Stroke(linewidth=2, foreground='black'), patheffects.Normal()])

            p2 = patches.Rectangle(
                (0.0, methylkeyori + 0.325),
                1.0,
                self.J_bCHG * 0.33,
                clip_on=True,
                lw=0,
                fill=False
            )
            ax2.add_patch(p2)
            chgc.set_clip_path(p2)
            chgh.set_clip_path(p2)
            chgg.set_clip_path(p2)

            # 绘制 CHH 图例
            CHHC = ax2.text(0.2, methylkeyori, 'C', color='white', ha='center', fontproperties=fontKey,
                            transform=transforms.offset_copy(
                                ax2.transData,
                                fig=self.fig,
                                x=0.0,
                                y=1.5,
                                units='points'
                            ))
            CHHH1 = ax2.text(0.5, methylkeyori, 'H', color='white', ha='center', fontproperties=fontKey,
                             transform=transforms.offset_copy(
                                 ax2.transData,
                                 fig=self.fig,
                                 x=0.0,
                                 y=1.5,
                                 units='points'
                             ))
            CHHH2 = ax2.text(0.8, methylkeyori, 'H', color='white', ha='center', fontproperties=fontKey,
                             transform=transforms.offset_copy(
                                 ax2.transData,
                                 fig=self.fig,
                                 x=0.0,
                                 y=1.5,
                                 units='points'
                             ))

            chhc = ax2.text(0.2, methylkeyori, 'C', color='black', ha='center', fontproperties=fontKey,
                            clip_on=True,
                            transform=transforms.offset_copy(
                                ax2.transData,
                                fig=self.fig,
                                x=0.0,
                                y=1.5,
                                units='points'
                            ))
            chhh1 = ax2.text(0.5, methylkeyori, 'H', color='black', ha='center', fontproperties=fontKey,
                             clip_on=True,
                             transform=transforms.offset_copy(
                                 ax2.transData,
                                 fig=self.fig,
                                 x=0.0,
                                 y=1.5,
                                 units='points'
                             ))
            chhh2 = ax2.text(0.8, methylkeyori, 'H', color='black', ha='center', fontproperties=fontKey,
                             clip_on=True,
                             transform=transforms.offset_copy(
                                 ax2.transData,
                                 fig=self.fig,
                                 x=0.0,
                                 y=1.5,
                                 units='points'
                             ))

            for txt in [CHHC, CHHH1, CHHH2, chhc, chhh1, chhh2]:
                txt.set_path_effects([patheffects.Stroke(linewidth=2, foreground='black'), patheffects.Normal()])

            p3 = patches.Rectangle(
                (0.0, methylkeyori),
                1.0,
                self.J_bCHH * 0.33,
                clip_on=True,
                lw=0,
                fill=False
            )
            ax2.add_patch(p3)
            chhc.set_clip_path(p3)
            chhh1.set_clip_path(p3)
            chhh2.set_clip_path(p3)

            # 绘制碱基频率图例
            freqkeyori = ax2top - 3
            bot = freqkeyori + 0.0
            for nt in sorted(['A', 'C', 'G', 'T'], reverse=True):
                freq = round(self.bgpps[nt], 4)
                ax2.bar(0.5, freq, bottom=bot, width=1.2, color=COLOR_SCHEME[nt], alpha=0.0)
                fontKey2 = FontProperties()
                if nt in ['C', 'G']:
                    fontKey2.set_size(60)
                    fontKey2.set_weight('bold')
                    fontKey2.set_family('monospace')
                    txt1 = ax2.text(
                        0.5,
                        bot,
                        nt,
                        ha='center',
                        va='baseline',
                        fontproperties=fontKey2,
                        transform=transforms.offset_copy(
                            ax2.transData,
                            fig=self.fig,
                            x=-5,
                            y=2.0 * self.bgpps[nt],
                            units='points'
                        ),
                        color=COLOR_SCHEME[nt],
                        zorder=1,
                    )
                    txt2 = ax2.text(
                        0.5,
                        bot,
                        nt,
                        ha='center',
                        va='baseline',
                        fontproperties=fontKey2,
                        transform=transforms.offset_copy(
                            ax2.transData,
                            fig=self.fig,
                            x=-5,
                            y=2.0 * self.bgpps[nt],
                            units='points'
                        ),
                        color='black',
                        clip_on=True,
                        zorder=2,
                    )
                    for txt in [txt1, txt2]:
                        txt.set_path_effects([Scale(1.2, 2.2 * self.bgpps[nt])])

                    bgMlevel = self.bgpps['CpG'] * self.J_bCG + self.bgpps['CHG'] * self.J_bCHG + self.bgpps['CHH'] * self.J_bCHH
                    pb = patches.Rectangle(
                        (0.0, bot),
                        1.0,
                        bgMlevel * 2.2 * self.bgpps[nt],
                        clip_on=True,
                        lw=0,
                        fill=False
                    )

                    ax2.add_patch(pb)
                    txt2.set_clip_path(pb)

                    ax2.plot(
                        (0.0, 1.0),
                        (bot + bgMlevel * 2.2 * self.bgpps[nt], bot + bgMlevel * 2.2 * self.bgpps[nt]),
                        linestyle='--',
                        linewidth=1,
                        color='lightgrey',
                    )
                else:
                    fontKey2.set_size(62)
                    fontKey2.set_weight('bold')
                    fontKey2.set_family('monospace')
                    txt = ax2.text(
                        0.5,
                        bot,
                        nt,
                        ha='center',
                        va='baseline',
                        fontproperties=fontKey2,
                        color=COLOR_SCHEME[nt]
                    )
                    txt.set_path_effects([Scale(1.0, 2.2 * self.bgpps[nt])])
                bot += freq * 2

            p3 = patches.Rectangle(
                (-0.1, freqkeyori - 0.05),
                1.2,
                3.1,
                clip_on=False,
                lw=1.5,
                fill=False
            )
            ax2.add_patch(p3)
            ax2.set_axis_off()
            ax2.text(
                0.5,
                -0.5,
                'MethylSeqLogo',
                ha='center',
                va='bottom',
                fontsize=10,
                weight='bold',
            )

            # 绘制甲基化信息熵在 ax3
            ax3height = 2 if Height == 3.0 else 5
            ax3.hlines(0, 0, self.plotlen - 1, color='blue', linewidth=0.5)
            ax3.set_ylim(-(ax3height / 2), ax3height / 2)
            ax3_y_axis = [-(int(ax3height / 2)), 0, int(ax3height / 2)]
            ax3.set_yticks(ax3_y_axis)
            ax3.set_yticklabels([abs(i) for i in ax3_y_axis])

            for i, cent in self.entropys.loc[:, self.entropys.columns.isin(
                    ['CpG_p', 'CpG_m', 'CHG_p', 'CHG_m', 'CHH_p', 'CHH_m'])].iterrows():
                Bottom_p = 0
                Bottom_m = 0
                cg_p, cg_m, chg_p, chg_m, chh_p, chh_m = cent

                ax3.bar(
                    x_axis[i],
                    cg_p,
                    bottom=Bottom_p,
                    width=0.8,
                    color='blue',
                    alpha=0.8,
                )
                Bottom_p += cg_p

                ax3.bar(
                    x_axis[i],
                    -cg_m,
                    bottom=Bottom_m,
                    width=0.8,
                    color='blue',
                    alpha=0.4,
                )
                Bottom_m -= cg_m

                ax3.bar(
                    x_axis[i],
                    chg_p,
                    bottom=Bottom_p,
                    width=0.8,
                    color='lime',
                    alpha=0.8,
                )
                Bottom_p += chg_p

                ax3.bar(
                    x_axis[i],
                    -chg_m,
                    bottom=Bottom_m,
                    width=0.8,
                    color='lime',
                    alpha=0.4,
                )
                Bottom_m -= chg_m

                ax3.bar(
                    x_axis[i],
                    chh_p,
                    bottom=Bottom_p,
                    width=0.8,
                    color='red',
                    alpha=0.8,
                )
                Bottom_p += chh_p

                ax3.bar(
                    x_axis[i],
                    -chh_m,
                    bottom=Bottom_m,
                    width=0.8,
                    color='red',
                    alpha=0.4,
                )
                Bottom_m += chh_m

            font = FontProperties(size=12, weight='bold', family='monospace')

            ax3.text(
                -0.5,
                max(ax3_y_axis),
                r'+ strand',
                color='black',
                ha='left',
                va='center',
                fontproperties=font,
            )
            ax3.text(
                -0.5,
                min(ax3_y_axis),
                r'- strand',
                color='black',
                ha='left',
                va='center',
                fontproperties=font,
            )

            ax3.spines['top'].set_linewidth(0)
            ax3.spines['right'].set_linewidth(0)
            ax3.spines['left'].set_linewidth(3)
            ax3.spines['bottom'].set_linewidth(0)
            sns.despine(ax=ax3, offset=2, trim=True)
            for label in ax3.yaxis.get_ticklabels():
                label.set_fontweight('bold')
            ax3.tick_params(
                direction='out',
                length=4,
                width=3,
                labelsize=18,
                top=False,
                right=False,
                bottom=False
            )
            ax3.axes.get_xaxis().set_visible(False)

            # 绘制甲基化信息熵图例在 ax4
            ax4height = 2 if ax3height == 2 else 2.5
            ax4.set_ylim(-(ax4height / 2), ax4height / 2)
            ax4_y_axis = [-(ax4height / 2), 0, ax4height / 2]
            ax4.set_yticks(ax4_y_axis)
            ax4.set_yticklabels([int(ax4height / 2), abs(0), abs(ax4height / 2)])

            ax4ori = 0 if ax4height == 2 else 0.25

            # 绘制 CHH、CHG、CG 图例
            p1 = patches.Rectangle(
                (0.0, ax4ori + 0.1),
                0.2,
                0.2,
                clip_on=False,
                lw=0,
                fill=True,
                color='red',
                alpha=0.8,
            )
            p2 = patches.Rectangle(
                (0.0, ax4ori + 0.35),
                0.2,
                0.2,
                clip_on=False,
                lw=0,
                fill=True,
                color='lime',
                alpha=0.8,
            )
            p3 = patches.Rectangle(
                (0.0, ax4ori + 0.60),
                0.2,
                0.2,
                clip_on=False,
                lw=0,
                fill=True,
                color='blue',
                alpha=0.8,
            )
            for p in [p1, p2, p3]:
                ax4.add_patch(p)

            ax4.text(
                0.35,
                ax4ori + 0.18,
                r'CHH',
                color='black',
                ha='left',
                va='center',
                fontproperties=font,
            )
            ax4.text(
                0.35,
                ax4ori + 0.43,
                r'CHG',
                color='black',
                ha='left',
                va='center',
                fontproperties=font,
            )
            ax4.text(
                0.35,
                ax4ori + 0.68,
                r'CG',
                color='black',
                ha='left',
                va='center',
                fontproperties=font,
            )

            frame = patches.Rectangle(
                (-0.1, ax4ori - 0.04),
                1.4,
                0.95,
                clip_on=False,
                lw=1.5,
                fill=False,
            )
            ax4.add_patch(frame)
            ax4.set_title('Entropy from:', fontsize=10, weight='bold')
            ax4.set_axis_off()

            # 绘制二聚体信息熵在 ax5
            ax5.set_ylim(-1, 3)
            ax5.hlines(0, 0, self.plotlen - 1, color='red', linewidth=0.3)

            ax5.bar([i + 0.5 for i in list(range(self.plotlen - 1))], self.dientropys, linewidth=2, color='red', alpha=0.6)
            ax5.hlines(self.bg_dientropys_max, 0, self.plotlen - 1, color='grey', linewidth=0.5, linestyle='--')
            ax5.hlines(self.bg_dientropys_min, 0, self.plotlen - 1, color='grey', linewidth=0.5, linestyle='--')

            ax5.set_yticks([-1, 0, 1, 2, 3])
            ax5.set_yticklabels([abs(i) for i in [-1, 0, 1, 2, 3]])
            ax5.spines['top'].set_linewidth(0)
            ax5.spines['right'].set_linewidth(0)
            ax5.spines['left'].set_linewidth(3)
            ax5.spines['bottom'].set_linewidth(0)

            sns.despine(ax=ax5, offset=2, trim=True)
            ax5.tick_params(
                direction='out',
                length=4,
                width=3,
                labelsize=18,
                top=False,
                right=False,
                bottom=False
            )
            for label in ax5.yaxis.get_ticklabels():
                label.set_fontweight('bold')
            ax5.axes.get_xaxis().set_visible(False)

            ax5.text(
                -0.5,
                3,
                r'+ dimer depleted',
                color='black',
                ha='left',
                va='center',
                fontproperties=font,
            )
            ax5.text(
                -0.5,
                -1,
                r'- dimer enriched',
                color='black',
                ha='left',
                va='center',
                fontproperties=font,
            )

            ax6.set_axis_off()
            ax6.text(
                0.0,
                self.bg_dientropys_max,
                r'max = ' + str(round(self.bg_dientropys_max, 2)),
                color='grey',
                ha='center',
                va='center',
                alpha=0.6,
                fontsize=10,
            )
            ax6.text(
                0.0,
                self.bg_dientropys_min,
                r'min = ' + str(round(self.bg_dientropys_min, 2)),
                color='grey',
                ha='center',
                va='center',
                alpha=0.6,
                fontsize=10,
            )

            # 设置 ax1 的样式
            ax1.set_yticks(y_axis)
            ax1.set_yticklabels([str(int(i)) for i in y_axis])
            ax1.spines['top'].set_linewidth(0)
            ax1.spines['right'].set_linewidth(0)
            ax1.spines['left'].set_linewidth(3)
            ax1.spines['bottom'].set_linewidth(3)

            sns.despine(ax=ax1, offset=2, trim=True)
            for label in ax1.xaxis.get_ticklabels():
                label.set_fontweight('bold')
            for label in ax1.yaxis.get_ticklabels():
                label.set_fontweight('bold')
            ax1.tick_params(
                direction='out',
                length=4,
                width=3,
                labelsize=18,
                top=False,
                right=False
            )

            ax1.set_ylabel('Bits', fontsize=18, weight='bold')
            ax1.set_title(self.TF, fontsize=16, weight='bold')

            return ax1, ax2, ax3, ax4, ax5

        else:
            ax1 = self.fig.add_subplot(111)

            x_axis = list(range(self.plotlen))
            ax1.set_xticks(x_axis)
            ax1.set_xticklabels([str(i + 1) for i in x_axis])

            y_axis = range(int(Height) + 1)

            font_T = FontProperties(size=70, weight='bold', family='monospace')
            font_C = FontProperties(size=68, weight='bold', family='monospace')

            bbox_props = dict(boxstyle="square, pad=0.0", fill=False, lw=0.0, alpha=0.5)

            base_heights = zip(self.base_heights, self.Cmethyls, self.Gmethyls)
            for i, (bases, _, _) in enumerate(base_heights):
                Bottom = 0
                width = 1

                for (base, score) in bases:
                    ax1.bar(
                        x_axis[i],
                        score,
                        bottom=Bottom,
                        width=width,
                        color=COLOR_SCHEME[base],
                        alpha=0.0
                    )

                    if base in ['C', 'G']:
                        yshift = 1.0
                        txt1 = ax1.text(
                            x_axis[i],
                            Bottom,
                            base,
                            fontproperties=font_C,
                            transform=transforms.offset_copy(
                                ax1.transData,
                                fig=self.fig,
                                x=-6,
                                y=yshift * score,
                                units='points'
                            ),
                            ha='center',
                            va='baseline',
                            color=COLOR_SCHEME[base],
                            bbox=bbox_props
                        )
                        txt1.set_path_effects([Scale(1.2, score)])
                    else:
                        txt1 = ax1.text(
                            x_axis[i],
                            Bottom,
                            base,
                            fontproperties=font_T,
                            ha='center',
                            va='baseline',
                            color=COLOR_SCHEME[base],
                            bbox=bbox_props
                        )
                        txt1.set_path_effects([Scale(1.0, score)])

                    Bottom += score

            ax1.set_yticks(y_axis)
            ax1.set_yticklabels([str(i) for i in y_axis])
            ax1.spines['top'].set_linewidth(0)
            ax1.spines['right'].set_linewidth(0)
            ax1.spines['left'].set_linewidth(3)
            ax1.spines['bottom'].set_linewidth(3)

            sns.despine(ax=ax1, offset=2, trim=True)
            for label in ax1.xaxis.get_ticklabels():
                label.set_fontweight('bold')
            for label in ax1.yaxis.get_ticklabels():
                label.set_fontweight('bold')
            ax1.tick_params(
                direction='out',
                length=4,
                width=3,
                labelsize=16,
                top=False,
                right=False
            )

            ax1.set_ylabel('Bits', fontsize=16, weight='bold')
            ax1.set_title(f"{self.TF} $\\cdot$ {self.celltype}", fontsize=18, weight='bold')
            ax1.text(
                -1.1,
                -0.5,
                'MethylSeqLogo',
                ha='center',
                va='bottom',
                fontsize=12,
                weight='bold',
            )
            return ax1,


class riverLake:
	def __init__(self, fig, tissue, fourletterppm, dippm, Cmethyls, Gmethyls, bgpps, J_bCG, J_bCHG, J_bCHH, Freqs_, mode, plotlen, TF,motiflen):
		self.fig = fig
		self.tissue = tissue
		self.fourletterppm = fourletterppm
		self.dippm = dippm
		self.Cmethyls = Cmethyls
		self.Gmethyls = Gmethyls
		self.bgpps = bgpps
		self.J_bCG = J_bCG
		self.J_bCHG = J_bCHG
		self.J_bCHH = J_bCHH
		self.Freqs_ = Freqs_
		self.mode = mode
		self.plotlen = plotlen
		self.TF = TF
		self.motiflen = motiflen
		
		

	def plotRiverLake(self):
	
		""" Plot methylriverlake """

		# remove axes content
		self.fig.clf()
		self.fig.tight_layout()

		Height = 3.0

		if (self.mode == 'Methyl'):
			gs = GridSpec(
							1, 
							2, 
							width_ratios=[self.motiflen-1, 1], 
							height_ratios = [1], 
							wspace = 0.0, 
							hspace = 0.0,
							)
			ax1 = plt.subplot(gs[0])
			ax2 = plt.subplot(gs[1])

			x_axis = list(range(self.plotlen))
			ax1.set_xticks(x_axis)
			ax1.set_xticklabels([str(i+1) for i in x_axis])
			ax1.set_ylim(-0.2, 1.2)
			
			bbox_props = dict(boxstyle = "square, pad = 0", fill = 0.0, lw = 0.0, alpha = 0.5)

			# convert pandas dataframe to numpy array with decimal = 1
			mfourletterppm = np.array(self.fourletterppm)
			mfourletterppm = np.round(mfourletterppm, 1)

			base_heights = zip(mfourletterppm, self.Cmethyls, self.Gmethyls)
			for i, j in enumerate(base_heights):

				# # xshift = 0
				# print i
				# print j
				yshift = 4
				Bottom = 0
				width = 1.0
				yindex = [0.0, 0.25, 0.5, 0.75]
				for (base, score) in zip(['A', 'C', 'G', 'T'], j[0]):
					ax1.bar(
							x_axis[i], 
							0.25,
							bottom = Bottom,
							width = width,
							color = COLOR_SCHEME[base],
							alpha = 0.0
							)

					if score > 0.00:
			
						oval = patches.Ellipse(
												xy = (x_axis[i], Bottom),
												width = 0.7 * score,
												height = 0.4 * score,
												color = 'lightgrey',
												alpha = 1.0,
												fill = True,
												zorder = 4
												)
						ax1.add_patch(oval)
						
						if i < (self.motiflen - 1):
							for k, (base1, score1) in enumerate(zip(['A', 'C', 'G', 'T'], mfourletterppm[i + 1])):
							
								if score1 > 0.00:
			
									index = base + base1
									ax1.plot([i, i+1], [Bottom, yindex[k]], lw = 25.0 * self.dippm[i][index], color = 'lightgrey', zorder = 2)

									if self.dippm[i][index] >= (score*score1):
										scale1 = self.dippm[i][index]/(score*score1)
										aline = ax1.plot([i, i+1], [Bottom, yindex[k]], lw = 1*scale1, color = 'navy', zorder = 3)
									else:
										pass

						# set font
						fontC = FontProperties()
						fontC.set_size(48 * score) # scale to probability
						fontC.set_weight('bold')
						fontC.set_family('monospace')

						txt1 = ax1.text(
										x_axis[i],
										Bottom,
										base,
										fontproperties = fontC,
										transform = transforms.offset_copy(
																			ax1.transData, 
																			fig = self.fig,
																			x = 0,
																			y = -yshift*score, 
																			units = 'points'
																			),
										ha = 'center',
										va = 'center',
										color = COLOR_SCHEME[base],
										bbox = bbox_props,
										zorder = 4,
										)
						if base in ['C', 'G']:
							txt2 = ax1.text(
											x_axis[i],
											Bottom,
											base,
											fontproperties = fontC,
											transform = transforms.offset_copy(
																				ax1.transData, 
																				fig = self.fig,
																				x = 0,
																				y = -yshift*score, 
																				units = 'points'
																				),
											ha = 'center',
											va = 'center',
											color = 'black',
											clip_on = True,
											bbox = bbox_props,
											zorder = 4,
											)
							txt3 = ax1.text(
											x_axis[i],
											Bottom,
											base,
											fontproperties = fontC,
											transform = transforms.offset_copy(
																				ax1.transData, 
																				fig = self.fig,
																				x = 0,
																				y = -yshift*score, 
																				units = 'points'
																				),
											ha = 'center',
											va = 'top',
											color = 'black',
											clip_on = True,
											alpha = 0.0,
											bbox = bbox_props,
											zorder = 4,
											)
							txt4 = ax1.text(
											x_axis[i],
											Bottom,
											base,
											fontproperties = fontC,
											transform = transforms.offset_copy(
																				ax1.transData, 
																				fig = self.fig,
																				x = 0,
																				y = -yshift*score, 
																				units = 'points'
																				),
											ha = 'center',
											va = 'baseline',
											color = 'black',
											clip_on = True,
											alpha = 0.0,
											bbox = bbox_props,
											zorder = 4,
											)

							self.fig.canvas.draw()
							bound1 = txt3.get_window_extent()
							bound2 = txt4.get_window_extent()

							inv = ax1.transData.inverted()
							a1 = inv.transform((bound1.xmin,bound1.ymin))
							a2 = inv.transform((bound1.xmax,bound1.ymax))
							a3 = inv.transform((bound2.xmin,bound2.ymin))
							a4 = inv.transform((bound2.xmax,bound2.ymax))
							
							letterheight = a4[1] - a2[1]
							
							left = x_axis[i] - 0.5
							halfletter = letterheight * 0.50
							
							bot = Bottom - halfletter

							bgMlevel = None
							foreMlevel = None
							if base == 'C':
								bgMlevel = self.Freqs_.iloc[i]['CpG_p'] * self.J_bCG + self.Freqs_.iloc[i]['CHG_p'] * self.J_bCHG + self.Freqs_.iloc[i]['CHH_p'] * self.J_bCHH
								foreMlevel = self.Freqs_.iloc[i]['CpG_p'] * j[1][0] + self.Freqs_.iloc[i]['CHG_p'] * j[1][1] + self.Freqs_.iloc[i]['CHH_p'] * j[1][2]
								
								# methylation = letterheight * j[1]
							elif base == 'G':
								bgMlevel = self.Freqs_.iloc[i]['CpG_m'] * self.J_bCG + self.Freqs_.iloc[i]['CHG_m'] * self.J_bCHG + self.Freqs_.iloc[i]['CHH_m'] * self.J_bCHH
								foreMlevel = self.Freqs_.iloc[i]['CpG_m'] * j[2][0] + self.Freqs_.iloc[i]['CHG_m'] * j[2][1] + self.Freqs_.iloc[i]['CHH_m'] * j[2][2]
								# methylation = letterheight * j[2]
							print (i, bot, base, letterheight, foreMlevel, bgMlevel)

							p = patches.Rectangle(
			    									(left, bot), 
			    									width,
			    									letterheight * foreMlevel, # height of the patch
			    									clip_on = True,
			    									lw = 0,
			    									fill = False,
			    									)					
							ax1.add_patch(p)
							txt2.set_clip_path(p)
							if letterheight >= 0.1:						
								ax1.plot(
											(i-0.3, i+0.3),
											(bot + letterheight * bgMlevel, bot + letterheight * bgMlevel),
											linestyle = '--',
											color = 'grey',
											linewidth = 1.0,
											zorder = 5, 
										)
						
					Bottom += 0.25

				# add key to logo
				fontKey = FontProperties()
				fontKey.set_size(22)
				fontKey.set_weight('bold')
				fontKey.set_family('monospace')
				
				ax2top = int(Height)
				methylkeyori = ax2top - 1
				
				# ax2.set_xlim(0.0, 1.0*(1+motiflen)/10)
				ax2.set_xlim(0.0, 1.2)
				ax2.set_ylim(0.0, ax2top)

				# y_axis = range(4)
				# ax2.set_yticks(y_axis)
				# ax2.set_yticklabels([str(i) for i in y_axis])
				# params = {'mathtext.default': 'regular' }  
				# plt.rcParams.update(params)			

				CGC = ax2.text(0.2, methylkeyori + 0.65, 'C', color = 'white', ha = 'center', fontproperties = fontKey,
							   transform = transforms.offset_copy(
								           							ax2.transData, 
								           							fig = self.fig,
								           							x = 0.0,
								           							y = 1.5, 
								           							units = 'points'
								           										),)
				CGG = ax2.text(0.5, methylkeyori + 0.65, 'G', color = 'white', ha = 'center', fontproperties = fontKey,
							   transform = transforms.offset_copy(
								           							ax2.transData, 
								           							fig = self.fig,
								           							x = 0.0,
								           							y = 1.5, 
								           							units = 'points'
								           										),)
				cgc = ax2.text(0.2, methylkeyori + 0.65, 'C', color = 'black', ha = 'center', fontproperties = fontKey, clip_on = True,
							   transform = transforms.offset_copy(
								           							ax2.transData, 
								           							fig = self.fig,
								           							x = 0.0,
								           							y = 1.5, 
								           							units = 'points'
								           										),)
				cgg = ax2.text(0.5, methylkeyori + 0.65, 'G', color = 'black', ha = 'center', fontproperties = fontKey, clip_on = True,
							   transform = transforms.offset_copy(
								           							ax2.transData, 
								           							fig = self.fig,
								           							x = 0.0,
								           							y = 1.5, 
								           							units = 'points'
								           										),)

				for txt in [CGC, CGG, cgc, cgg]:
					txt.set_path_effects([matplotlib.patheffects.Stroke(linewidth=2, foreground='black'), matplotlib.patheffects.Normal()])
					

				p1 = patches.Rectangle(
										(0.0, methylkeyori + 0.65), 
										1.0, 
										self.J_bCG*0.33, # height of the patch
										clip_on = True,
										lw = 0,
										fill = False
										)					
				ax2.add_patch(p1)
				cgc.set_clip_path(p1)
				cgg.set_clip_path(p1)


				CHGC = ax2.text(0.2, methylkeyori + 0.325, 'C', color = 'white', ha = 'center', fontproperties = fontKey, 
							   transform = transforms.offset_copy(
								           							ax2.transData, 
								           							fig = self.fig,
								           							x = 0.0,
								           							y = 1.5, 
								           							units = 'points'
								           										),)
				CHGH = ax2.text(0.5, methylkeyori + 0.325, 'H', color = 'white', ha = 'center', fontproperties = fontKey,
							   transform = transforms.offset_copy(
								           							ax2.transData, 
								           							fig = self.fig,
								           							x = 0.0,
								           							y = 1.5, 
								           							units = 'points'
								           										),)
				CHGG = ax2.text(0.8, methylkeyori + 0.325, 'G', color = 'white', ha = 'center', fontproperties = fontKey,
							   transform = transforms.offset_copy(
								           							ax2.transData, 
								           							fig = self.fig,
								           							x = 0.0,
								           							y = 1.5, 
								           							units = 'points'
								           										),)

				chgc = ax2.text(0.2, methylkeyori + 0.325, 'C', color = 'black', ha = 'center', fontproperties = fontKey, clip_on = True,
							   transform = transforms.offset_copy(
								           							ax2.transData, 
								           							fig = self.fig,
								           							x = 0.0,
								           							y = 1.5, 
								           							units = 'points'
								           										),)
				chgh = ax2.text(0.5, methylkeyori + 0.325, 'H', color = 'black', ha = 'center', fontproperties = fontKey, clip_on = True,
							   transform = transforms.offset_copy(
								           							ax2.transData, 
								           							fig = self.fig,
								           							x = 0.0,
								           							y = 1.5, 
								           							units = 'points'
								           										),)
				chgg = ax2.text(0.8, methylkeyori + 0.325, 'G', color = 'black', ha = 'center', fontproperties = fontKey, clip_on = True,
							   transform = transforms.offset_copy(
								           							ax2.transData, 
								           							fig = self.fig,
								           							x = 0.0,
								           							y = 1.5, 
								           							units = 'points'
								           										),)

				for txt in [CHGC, CHGH, CHGG, chgc, chgh, chgg]:
					txt.set_path_effects([matplotlib.patheffects.Stroke(linewidth=2, foreground='black'), matplotlib.patheffects.Normal()])

				p2 = patches.Rectangle(
										(0.0, methylkeyori + 0.325), 
										1.0, 
										self.J_bCHG*0.33, # height of the patch
										clip_on = True,
										lw = 0,
										fill = False
										)					
				ax2.add_patch(p2)
				chgc.set_clip_path(p2)
				chgh.set_clip_path(p2)
				chgg.set_clip_path(p2)

				CHHC = ax2.text(0.2, methylkeyori, 'C', color = 'white', ha = 'center', fontproperties = fontKey, 
							   transform = transforms.offset_copy(
								           							ax2.transData, 
								           							fig = self.fig,
								           							x = 0.0,
								           							y = 1.5, 
								           							units = 'points'
								           										),)
				CHHH1 = ax2.text(0.5, methylkeyori, 'H', color = 'white', ha = 'center', fontproperties = fontKey,
							   transform = transforms.offset_copy(
								           							ax2.transData, 
								           							fig = self.fig,
								           							x = 0.0,
								           							y = 1.5, 
								           							units = 'points'
								           										),)
				CHHH2 = ax2.text(0.8, methylkeyori, 'H', color = 'white', ha = 'center', fontproperties = fontKey,
							   transform = transforms.offset_copy(
								           							ax2.transData, 
								           							fig = self.fig,
								           							x = 0.0,
								           							y = 1.5, 
								           							units = 'points'
								           										),)

				chhc = ax2.text(0.2, methylkeyori, 'C', color = 'black', ha = 'center', fontproperties = fontKey, clip_on = True,
							   transform = transforms.offset_copy(
								           							ax2.transData, 
								           							fig = self.fig,
								           							x = 0.0,
								           							y = 1.5, 
								           							units = 'points'
								           										),)
				chhh1 = ax2.text(0.5, methylkeyori, 'H', color = 'black', ha = 'center', fontproperties = fontKey, clip_on = True,
							   transform = transforms.offset_copy(
								           							ax2.transData, 
								           							fig = self.fig,
								           							x = 0.0,
								           							y = 1.5, 
								           							units = 'points'
								           										),)
				chhh2 = ax2.text(0.8, methylkeyori, 'H', color = 'black', ha = 'center', fontproperties = fontKey, clip_on = True,
							   transform = transforms.offset_copy(
								           							ax2.transData, 
								           							fig = self.fig,
								           							x = 0.0,
								           							y = 1.5, 
								           							units = 'points'
								           										),)

				for txt in [CHHC, CHHH1, CHHH2, chhc, chhh1, chhh2]:
					# txt.set_path_effects([Scale(1.0, self.bgpps[nt])])
					txt.set_path_effects([matplotlib.patheffects.Stroke(linewidth=2, foreground='black'), matplotlib.patheffects.Normal()])

				p3 = patches.Rectangle(
										(0.0, methylkeyori), 
										1.0, 
										self.J_bCHH*0.33, # height of the patch
										clip_on = True,
										lw = 0,
										fill = False
										)					
				ax2.add_patch(p3)
				chhc.set_clip_path(p3)
				chhh1.set_clip_path(p3)
				chhh2.set_clip_path(p3)


				# frequency key
				freqkeyori = ax2top - 3.0

				bot = freqkeyori + 0.2
				for nt in sorted(['A', 'C', 'G', 'T'], reverse = True):
					freq = round(self.bgpps[nt], 4)
					bar = ax2.bar(0.5, freq, bottom = bot, width = 1.0, color = COLOR_SCHEME[nt], alpha = 0.0)
					oval = patches.Ellipse(
												xy = (0.5, bot),
												# radius = 1.0 * freq * 0.5,
												width = 1.0 * freq,
												height = 0.9 * freq,
												color = 'lightgrey',
												alpha = 1.0,
												fill = True,
												clip_on = False,
												zorder = 3,
												)
					ax2.add_patch(oval)

					fontKey2 = FontProperties()
					fontKey2.set_size(50 * self.bgpps[nt])
					fontKey2.set_weight('bold')
					fontKey2.set_family('monospace')
					
					txt = ax2.text(
										0.5, 
										bot, 
										nt, 
										ha = 'center', 
										va = 'center', 
										fontproperties = fontKey2,
										color = COLOR_SCHEME[nt],
										zorder = 4,
										)
					if nt in ['C', 'G']:
						txt2 = ax2.text(
											0.5, 
											bot, 
											nt, 
											ha = 'center', 
											va = 'center', 
											fontproperties = fontKey2,
											color = 'black',
											clip_on = True,
											zorder = 4,
											)
						txt3 = ax2.text(
											0.5, 
											bot, 
											nt, 
											ha = 'center', 
											va = 'top', 
											fontproperties = fontKey2,
											color = 'black',
											alpha = 0.0,
											)
						txt4 = ax2.text(
											0.5, 
											bot, 
											nt, 
											ha = 'center', 
											va = 'baseline', 
											fontproperties = fontKey2,
											color = 'black',
											alpha = 0.0,
											)
						self.fig.canvas.draw()
						bound1 = txt3.get_window_extent()
						bound2 = txt4.get_window_extent()

						inv = ax2.transData.inverted()
						a1 = inv.transform((bound1.xmin,bound1.ymin))
						a2 = inv.transform((bound1.xmax,bound1.ymax))
						a3 = inv.transform((bound2.xmin,bound2.ymin))
						a4 = inv.transform((bound2.xmax,bound2.ymax))
							
						letterheight = a4[1] - a2[1]
						bgMlevel = None
						bgMlevel = self.bgpps['CpG'] * self.J_bCG + self.bgpps['CHG'] * self.J_bCHG + self.bgpps['CHH'] * self.J_bCHH
						halfletterheight = letterheight/2.0
	
						pb = patches.Rectangle(
												(0.0, bot - halfletterheight),
												1.0,
												bgMlevel*letterheight,
												clip_on = True,
												lw = 0.0,
												fill = False,
												)
						ax2.add_patch(pb)
						txt2.set_clip_path(pb)
						ax2.plot(
									(0.3, 0.7),
									(bot - halfletterheight + bgMlevel*letterheight, bot - halfletterheight + bgMlevel*letterheight),
									linestyle = '--',
									linewidth = 0.5,
									color = 'grey',
									zorder = 5,
									)

					bot = bot + 0.25*2

				p31 = patches.Rectangle(
										(-0.1, freqkeyori-0.05), 
										1.2, 
										3.1, # height of the patch
										clip_on = False,
										lw = 1.2,
										fill = False
										)
				ax2.add_patch(p31)
				ax2.set_axis_off()
				ax2.text(	
							0.5, 
							-0.5, 
							'MethylSeqLogo', 
							ha = 'center', 
							va = 'bottom',
							fontsize = 10,
							weight = 'bold',
							)

			ax1.spines['top'].set_linewidth(0)
			ax1.spines['right'].set_linewidth(0)
			ax1.spines['left'].set_linewidth(0)
			ax1.spines['bottom'].set_linewidth(2)

			seaborn.despine(ax = ax1, offset = 2, trim = True)		
			for label in ax1.xaxis.get_ticklabels():
				label.set_fontweight('bold')

			ax1.get_yaxis().set_visible(False)
			ax1.tick_params(direction='out', length=4, width=2, labelsize= 'x-large', 
							top=False, right=False, left = False)
			# ax1.text(motiflen, -0.4, 'MethylSeqLogo', ha = 'center', va = 'bottom', fontsize = 10)
			# ax1.set_ylabel('Nucleotide', fontsize = 14, weight = 'bold')
			ax1.set_title(self.TF, fontsize = 16, weight = 'bold')

			return ax1, ax2,
		

		else:

			ax1 = self.fig.add_subplot(111)
			x_axis = list(range(self.plotlen))
			ax1.set_xticks(x_axis)
			ax1.set_xticklabels([str(i+1) for i in x_axis])
			ax1.set_ylim(-0.2, 1.2)
			
			# y_axis = range(int(Height) + 1)

			bbox_props = dict(boxstyle = "square, pad = 0", fill = 0.0, lw = 0.0, alpha = 0.5)

			# convert pandas dataframe to numpy array with decimal = 1
			mfourletterppm = np.array(self.fourletterppm)
			mfourletterppm = np.round(mfourletterppm, 1)

			# base_heights = zip(mfourletterppm, self.Cmethyls, self.Gmethyls)
			for i, j in enumerate(mfourletterppm):
				print (i, j)
				# xshift = 0
				yshift = 4
				Bottom = 0
				width = 1.0
				yindex = [0.0, 0.25, 0.5, 0.75]
				for (base, score) in zip(['A', 'C', 'G', 'T'], j):
					ax1.bar(
							x_axis[i], 
							0.25,
							bottom = Bottom,
							width = width,
							color = COLOR_SCHEME[base],
							alpha = 0.0
							)

					if score > 0.00:
			
						oval = patches.Ellipse(
												xy = (x_axis[i], Bottom),
												width = 0.7 * score,
												height = 0.4 * score,
												color = 'lightgrey',
												alpha = 1.0,
												fill = True,
												zorder = 3
												)
						ax1.add_patch(oval)
						
						if i < (self.motiflen - 1):
							for k, (base1, score1) in enumerate(zip(['A', 'C', 'G', 'T'], mfourletterppm[i + 1])):
							
								if score1 > 0.00:
			
									index = base + base1
									ax1.plot([i, i+1], [Bottom, yindex[k]], lw = 25.0 * self.dippm[i][index], color = 'lightgrey', zorder = 1)

									if self.dippm[i][index] >= (score*score1):
										scale1 = self.dippm[i][index]/(score*score1)
										aline = ax1.plot([i, i+1], [Bottom, yindex[k]], lw = 1*scale1, color = 'navy', zorder = 1)
									else:
										pass

						# set font
						fontC = FontProperties()
						fontC.set_size(48 * score) # scale to probability
						fontC.set_weight('bold')
						fontC.set_family('monospace')

						txt1 = ax1.text(
										x_axis[i],
										Bottom,
										base,
										fontproperties = fontC,
										transform = transforms.offset_copy(
																			ax1.transData, 
																			fig = self.fig,
																			x = 0,
																			y = -yshift*score, 
																			units = 'points'
																			),
										ha = 'center',
										va = 'center',
										color = COLOR_SCHEME[base],
										bbox = bbox_props
										)
					Bottom += 0.25
			ax1.spines['top'].set_linewidth(0)
			ax1.spines['right'].set_linewidth(0)
			ax1.spines['left'].set_linewidth(0)
			ax1.spines['bottom'].set_linewidth(2)

			seaborn.despine(ax = ax1, offset = 2, trim = True)		
			for label in ax1.xaxis.get_ticklabels():
				label.set_fontweight('bold')

			ax1.get_yaxis().set_visible(False)
			ax1.tick_params(
							direction='out', 
							length=4, 
							width=2, 
							labelsize= 16,#'x-large', 
							top=False, 
							right=False, 
							left = False)
			ax1.text(
						self.motiflen-1, 
						-0.6, 
						'MethylSeqLogo', 
						ha = 'center', 
						va = 'bottom',
						fontsize = 12,
						weight = 'bold',
						)
			# ax1.set_ylabel('Nucleotide', fontsize = 14, weight = 'bold')
			ax1.set_title(self.TF, fontsize = 18, weight = 'bold')


			return ax1,