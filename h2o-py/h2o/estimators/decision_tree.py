#!/usr/bin/env python
# -*- encoding: utf-8 -*-
#
# This file is auto-generated by h2o-3/h2o-bindings/bin/gen_python.py
# Copyright 2016 H2O.ai;  Apache License Version 2.0 (see LICENSE for details)
#

from h2o.estimators.estimator_base import H2OEstimator
from h2o.exceptions import H2OValueError
from h2o.frame import H2OFrame
from h2o.utils.typechecks import assert_is_type, Enum, numeric


class H2ODecisionTreeEstimator(H2OEstimator):
    """
    Decision Tree

    """

    algo = "dt"
    supervised_learning = True

    def __init__(self,
                 model_id=None,  # type: Optional[Union[None, str, H2OEstimator]]
                 training_frame=None,  # type: Optional[Union[None, str, H2OFrame]]
                 ignored_columns=None,  # type: Optional[List[str]]
                 ignore_const_cols=True,  # type: bool
                 categorical_encoding="auto",  # type: Literal["auto", "enum", "one_hot_internal", "one_hot_explicit", "binary", "eigen", "label_encoder", "sort_by_response", "enum_limited"]
                 response_column=None,  # type: Optional[str]
                 seed=-1,  # type: int
                 max_depth=20,  # type: int
                 min_rows=10,  # type: int
                 ):
        """
        :param model_id: Destination id for this model; auto-generated if not specified.
               Defaults to ``None``.
        :type model_id: Union[None, str, H2OEstimator], optional
        :param training_frame: Id of the training data frame.
               Defaults to ``None``.
        :type training_frame: Union[None, str, H2OFrame], optional
        :param ignored_columns: Names of columns to ignore for training.
               Defaults to ``None``.
        :type ignored_columns: List[str], optional
        :param ignore_const_cols: Ignore constant columns.
               Defaults to ``True``.
        :type ignore_const_cols: bool
        :param categorical_encoding: Encoding scheme for categorical features
               Defaults to ``"auto"``.
        :type categorical_encoding: Literal["auto", "enum", "one_hot_internal", "one_hot_explicit", "binary", "eigen", "label_encoder",
               "sort_by_response", "enum_limited"]
        :param response_column: Response variable column.
               Defaults to ``None``.
        :type response_column: str, optional
        :param seed: Seed for random numbers (affects sampling)
               Defaults to ``-1``.
        :type seed: int
        :param max_depth: Max depth of tree.
               Defaults to ``20``.
        :type max_depth: int
        :param min_rows: Fewest allowed (weighted) observations in a leaf.
               Defaults to ``10``.
        :type min_rows: int
        """
        super(H2ODecisionTreeEstimator, self).__init__()
        self._parms = {}
        self._id = self._parms['model_id'] = model_id
        self.training_frame = training_frame
        self.ignored_columns = ignored_columns
        self.ignore_const_cols = ignore_const_cols
        self.categorical_encoding = categorical_encoding
        self.response_column = response_column
        self.seed = seed
        self.max_depth = max_depth
        self.min_rows = min_rows

    @property
    def training_frame(self):
        """
        Id of the training data frame.

        Type: ``Union[None, str, H2OFrame]``.
        """
        return self._parms.get("training_frame")

    @training_frame.setter
    def training_frame(self, training_frame):
        self._parms["training_frame"] = H2OFrame._validate(training_frame, 'training_frame')

    @property
    def ignored_columns(self):
        """
        Names of columns to ignore for training.

        Type: ``List[str]``.
        """
        return self._parms.get("ignored_columns")

    @ignored_columns.setter
    def ignored_columns(self, ignored_columns):
        assert_is_type(ignored_columns, None, [str])
        self._parms["ignored_columns"] = ignored_columns

    @property
    def ignore_const_cols(self):
        """
        Ignore constant columns.

        Type: ``bool``, defaults to ``True``.
        """
        return self._parms.get("ignore_const_cols")

    @ignore_const_cols.setter
    def ignore_const_cols(self, ignore_const_cols):
        assert_is_type(ignore_const_cols, None, bool)
        self._parms["ignore_const_cols"] = ignore_const_cols

    @property
    def categorical_encoding(self):
        """
        Encoding scheme for categorical features

        Type: ``Literal["auto", "enum", "one_hot_internal", "one_hot_explicit", "binary", "eigen", "label_encoder",
        "sort_by_response", "enum_limited"]``, defaults to ``"auto"``.
        """
        return self._parms.get("categorical_encoding")

    @categorical_encoding.setter
    def categorical_encoding(self, categorical_encoding):
        assert_is_type(categorical_encoding, None, Enum("auto", "enum", "one_hot_internal", "one_hot_explicit", "binary", "eigen", "label_encoder", "sort_by_response", "enum_limited"))
        self._parms["categorical_encoding"] = categorical_encoding

    @property
    def response_column(self):
        """
        Response variable column.

        Type: ``str``.
        """
        return self._parms.get("response_column")

    @response_column.setter
    def response_column(self, response_column):
        assert_is_type(response_column, None, str)
        self._parms["response_column"] = response_column

    @property
    def seed(self):
        """
        Seed for random numbers (affects sampling)

        Type: ``int``, defaults to ``-1``.
        """
        return self._parms.get("seed")

    @seed.setter
    def seed(self, seed):
        assert_is_type(seed, None, int)
        self._parms["seed"] = seed

    @property
    def max_depth(self):
        """
        Max depth of tree.

        Type: ``int``, defaults to ``20``.
        """
        return self._parms.get("max_depth")

    @max_depth.setter
    def max_depth(self, max_depth):
        assert_is_type(max_depth, None, int)
        self._parms["max_depth"] = max_depth

    @property
    def min_rows(self):
        """
        Fewest allowed (weighted) observations in a leaf.

        Type: ``int``, defaults to ``10``.
        """
        return self._parms.get("min_rows")

    @min_rows.setter
    def min_rows(self, min_rows):
        assert_is_type(min_rows, None, int)
        self._parms["min_rows"] = min_rows


