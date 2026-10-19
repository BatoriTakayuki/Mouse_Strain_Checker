"""Conservative reference scores after explicitly grouping the two B6 strains.

This summary changes only the displayed candidate resolution. A group score is
its best exact model score, and is neither a probability nor evidence of F1,
mixture, or a particular unmeasured substrain.
"""


_STRAIN_GROUPS = {
    "C57BL_6J": "C57BL/6 (6J/6NJ)",
    "C57BL_6NJ": "C57BL/6 (6J/6NJ)",
}


def background_result(summary, rules, window=False):
    """Return a descriptive background group alongside the unchanged exact call.

    The two explicitly named B6 candidates share a display group only when both
    are present in this panel. The winning exact model must clear both original
    separation requirements against *every* model outside its group. Members
    within a group are not scored collectively.
"""
    row = dict(background_status="insufficient_sites", background_call="",
               background_candidate="", background_runner_up="", background_members="",
               background_diagnostic_sites=0, background_margin_per_site=None,
               background_fit_loss_per_site=None, background_reason="背景参考評価に必要な観測SNP数に達していません。")
    if not summary.n:
        row.update(background_status="no_observed_snps", background_reason="評価に使えるSNPが観測されていません。")
        return row

    strain_groups = _STRAIN_GROUPS if set(_STRAIN_GROUPS).issubset(summary.strains) else {}
    groups = {}
    model_groups = []
    for index, (a, b) in enumerate(summary.models):
        group = tuple(sorted({strain_groups.get(summary.strains[a], summary.strains[a]),
                              strain_groups.get(summary.strains[b], summary.strains[b])}))
        groups.setdefault(group, []).append(index)
        model_groups.append(group)

    top = max(range(len(summary.models)), key=summary.scores.__getitem__)
    best_group = model_groups[top]
    outsiders = [i for i, group in enumerate(model_groups) if group != best_group]
    loss = max(0.0, (summary.saturated - summary.scores[top]) / summary.n)
    row.update(background_candidate=" + ".join(best_group),
               background_members=" | ".join(summary.name(i) for i in groups[best_group]),
               background_fit_loss_per_site=loss)

    if outsiders:
        second = max(outsiders, key=summary.scores.__getitem__)
        a, b = summary.models[top]
        diagnostics = []
        for index in outsiders:
            c, d = summary.models[index]
            diagnostics.append(sum(n for genotypes, n in summary.patterns.items()
                                   if int(genotypes[a]) + int(genotypes[b])
                                   != int(genotypes[c]) + int(genotypes[d])))
        diagnostic = min(diagnostics)
        margin = (summary.scores[top] - summary.scores[second]) / summary.n
        row.update(background_runner_up=" + ".join(model_groups[second]),
                   background_diagnostic_sites=diagnostic, background_margin_per_site=margin)

    if summary.n < (rules.min_window_sites if window else rules.min_sites):
        return row
    if loss > rules.max_fit_loss:
        row.update(background_status="poor_fit", background_reason="最上位モデルでも観測SNPへの適合が基準を満たしていません。")
    elif not outsiders:
        row.update(background_status="ambiguous", background_reason="比較できる別の背景候補群がないため、背景を区別できません。")
    else:
        reasons = []
        if diagnostic < rules.min_diagnostic:
            reasons.append("観測SNP上で別の背景候補群と区別できません。" if diagnostic == 0
                           else f"別の背景候補群と区別できるSNPが{diagnostic}個で、基準{rules.min_diagnostic}個に達していません。")
        if margin < rules.min_margin:
            reasons.append("別の背景候補群とのスコア差が基準未満です。")
        if reasons:
            row.update(background_status="ambiguous", background_reason="".join(reasons))
        else:
            row.update(background_status="candidate_group_supported",
                       background_call=row["background_candidate"], background_reason="")
    return row
