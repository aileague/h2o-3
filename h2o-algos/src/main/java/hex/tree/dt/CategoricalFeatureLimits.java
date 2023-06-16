package hex.tree.dt;

import java.util.Arrays;

/**
 * Limits for one feature.
 */
public class CategoricalFeatureLimits extends AbstractFeatureLimits {
    public int[] _setOfCategories;
    public CategoricalFeatureLimits(final int[] setOfCategories) {
        _setOfCategories = Arrays.copyOf(setOfCategories, setOfCategories.length);
    }

    public void setNewSetOfCategories(final int[] newSetOfCategories) {
        _setOfCategories = Arrays.copyOf(newSetOfCategories, newSetOfCategories.length);
    }

    public CategoricalFeatureLimits clone() {
        return new CategoricalFeatureLimits(Arrays.copyOf(_setOfCategories, _setOfCategories.length));
    }

    @Override
    public double[] toDoubles() {
        return Arrays.stream(_setOfCategories).mapToDouble(c -> (double) c).toArray();
    }

    @Override
    public boolean equals(AbstractFeatureLimits other) {
        return Arrays.equals(_setOfCategories, ((CategoricalFeatureLimits) other)._setOfCategories);
    }
}
