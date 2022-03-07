package com.hartwig.hmftools.compar;

import static com.hartwig.hmftools.compar.CommonUtils.DATA_DELIM;

import java.util.StringJoiner;

public class MismatchOld
{
    public final Category Category;
    public final MismatchType MismatchType;
    public final boolean Reportable;
    public final String RefSource;
    public final String OtherSource;
    public final String ItemName;
    public final String Gene; // if applicable
    public final String DiffValue;

    public MismatchOld(
            final Category category, final MismatchType mismatchType, final String refSource, final String otherSource,
            final boolean reportable, final String itemName, final String gene, final String diffValue)
    {
        Category = category;
        MismatchType = mismatchType;
        Reportable = reportable;
        RefSource = refSource;
        OtherSource = otherSource;
        ItemName = itemName;
        Gene = gene;
        DiffValue = diffValue;
    }

    public static String header()
    {
        return "Category,MismatchType,Reportable,RefSource,OtherSource,Item,Gene,Diff";
    }

    // Category,Reportable,RefItem,RefSource,OtherSource,RefValues,OtherValues,RefDiffs,OtherDiffs

    public String toCsv()
    {
        StringJoiner sj = new StringJoiner(DATA_DELIM);
        sj.add(Category.toString());
        sj.add(MismatchType.toString());
        sj.add(String.valueOf(Reportable));
        sj.add(RefSource);
        sj.add(OtherSource);
        sj.add(ItemName);
        sj.add(Gene);
        sj.add(DiffValue);
        return sj.toString();
    }
}