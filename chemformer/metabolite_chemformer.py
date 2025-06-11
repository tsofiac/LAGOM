import logging
import pandas as pd
import torch
import numpy as np
from typing import Any, Dict, List, Optional, Tuple, Union

from aizynthmodels.chemformer.data.encoder import BatchEncoder
from aizynthmodels.chemformer.data.tokenizer import ReplaceTokensMasker, TokensMasker
from aizynthmodels.chemformer.data.base import BaseDataModule
from aizynthmodels.utils.tokenizer import SMILESAugmenter
from aizynthmodels.utils.type_utils import PredictionType, StrDict, TargetType
from aizynthmodels.utils.smiles import canonicalize_smiles
from aizynthmodels.utils.scores import BaseScore


class MulitpleReactionListDataModule(BaseDataModule):
    def __init__(
        self,
        augment_prob: float = 0.0,
        reverse: bool = False,
        masker: Optional[Union[TokensMasker, ReplaceTokensMasker]] = None,
        **kwargs,
    ):
        super().__init__(**kwargs)
        self._batch_augmenter = SMILESAugmenter(augment_prob=augment_prob)
        self._encoder = BatchEncoder(
            tokenizer=self.tokenizer, masker=masker, max_seq_len=self.max_seq_len
        )
        self.reverse = reverse

    def _collate(self, batch: List[StrDict], train: bool = True) -> StrDict:
        (
            encoder_ids,
            encoder_mask,
            decoder_ids,
            decoder_mask,
            input_smiles,
            target_smiles,
            possible_products,
        ) = self._transform_batch(batch, train)

        return {
            "encoder_input": encoder_ids,
            "encoder_pad_mask": encoder_mask,
            "decoder_input": decoder_ids[:-1, :],
            "decoder_pad_mask": decoder_mask[:-1, :],
            "target": decoder_ids.clone()[1:, :],
            "target_mask": decoder_mask.clone()[1:, :],
            "target_smiles": possible_products,  # target_smiles
            "input_smiles": input_smiles,
        }

    def _get_sequences(
        self, batch: List[StrDict], train: bool
    ) -> Tuple[List[str], List[str]]:
        reactants = [item["reactants"] for item in batch]
        products = [item["products"] for item in batch]

        if train and self._batch_augmenter.augment_prob > 0.0:
            reactants = self._batch_augmenter(reactants)
            products = self._batch_augmenter(products)

        return reactants, products

    def _load_all_data(self) -> None:
        with open(self.dataset_path, "r") as fileobj:
            lines = fileobj.read().splitlines()
        if ">>" in lines[0]:
            reactants, products = zip(*[line.split(">>") for line in lines])
        else:
            reactants = lines
            products = lines.copy()
        self._all_data = {"reactants": reactants, "products": products}

    def _transform_batch(
        self, batch: List[StrDict], train: bool
    ) -> Tuple[torch.Tensor, torch.Tensor, torch.Tensor, torch.Tensor, List[str]]:
        reactants_smiles, product_smiles, possible_products = self._get_sequences(
            batch, train
        )
        reactants_ids, reactants_mask = self._encoder(
            reactants_smiles,
            add_sep_token=self.add_sep_token and not self.reverse,
            mask=not self.reverse,
        )
        product_ids, product_mask = self._encoder(
            product_smiles,
            add_sep_token=self.add_sep_token and self.reverse,
            mask=self.reverse,
        )
        if not self.reverse:
            return (
                reactants_ids,
                reactants_mask,
                product_ids,
                product_mask,
                reactants_smiles,
                product_smiles,
                possible_products,
            )
        return (
            product_ids,
            product_mask,
            reactants_ids,
            reactants_mask,
            product_smiles,
            reactants_smiles,
            possible_products,
        )


class MetaboliteDataModule(MulitpleReactionListDataModule):
    """
    For using InPossibleProductsScore for metabolite prediction
    """

    def __init__(
        self,
        reactants: Optional[List[str]] = None,
        products: Optional[List[str]] = None,
        possible_products: Optional[List[str]] = None,
        augmentation_strategy: Optional[str] = None,
        **kwargs,
    ):
        super().__init__(**kwargs)

        self._augment_strategy = augmentation_strategy
        self._in_memory = False
        if reactants is not None and products is not None:
            self._in_memory = True
            logging.info("Using in-memory datamodule.")
            self._all_data = {
                "reactants": reactants,
                "products": products,
                "possible_products": possible_products,
            }

    def _get_sequences(
        self, batch: List[Dict[str, Any]], train: bool
    ) -> Tuple[List[str], List[str]]:
        reactants = [item["reactants"] for item in batch]
        products = [item["products"] for item in batch]
        possible_products = [item["possible_products"] for item in batch]

        if train:
            if self._augment_strategy == "reactants" or self._augment_strategy == "all":
                reactants = self._batch_augmenter(reactants)
            if self._augment_strategy == "products" or self._augment_strategy == "all":
                products = self._batch_augmenter(products)
        return reactants, products, possible_products

    def _load_all_data(self) -> None:
        if self._in_memory:
            return

        if self.dataset_path.endswith(".csv"):
            df = pd.read_csv(self.dataset_path, sep="\t").reset_index()
            self._all_data = {
                "reactants": df["reactants"].tolist(),
                "products": df["products"].tolist(),
                "possible_products": df["possible_products"].tolist(),
            }
            self._set_split_indices_from_dataframe(df)
        else:
            super()._load_all_data()


class Top1MetaboliteScore(BaseScore):
    scorer_name = "top_1_metabolite_score"

    def __init__(
        self,
        canonicalized: bool = False,
    ):
        """
        Uses the MetaboliteDataModule
        :param canonicalized: whether the sampled_smiles and target_smiles are
                been canonicalized.
        """
        super().__init__()
        self._canonicalized = canonicalized

    def _is_in_set(
        self, predictions: PredictionType, ground_truth: TargetType
    ) -> np.ndarray:
        ground_truth = [smiles.split(".") for smiles in ground_truth]

        if not self._canonicalized:
            predictions = [
                [canonicalize_smiles(smiles) for smiles in smiles_list]
                for smiles_list in predictions
            ]
            ground_truth = [
                [canonicalize_smiles(smiles) for smiles in smiles_list]
                for smiles_list in ground_truth
            ]

        is_in_set = [
            sampled_smi[0] in tgt_smi if len(tgt_smi) > 0 else False
            for sampled_smi, tgt_smi in zip(predictions, ground_truth)
        ]

        return is_in_set

    def _score_predictions(
        self, predictions: PredictionType, ground_truth: TargetType
    ) -> Dict[str, float]:
        is_in_set = self._is_in_set(predictions, ground_truth)

        score = np.mean(is_in_set)

        return {self.scorer_name: score}


class MetaboliteCoverageScore(BaseScore):
    scorer_name = "metabolite_coverage_score"

    def __init__(
        self,
        canonicalized: bool = False,
    ):
        """
        Uses the MetaboliteDataModule
        :param canonicalized: whether the sampled_smiles and target_smiles are
                been canonicalized.
        """
        super().__init__()
        self._canonicalized = canonicalized

    def _is_in_set(
        self, predictions: PredictionType, ground_truth: TargetType
    ) -> np.ndarray:
        ground_truth = [smiles.split(".") for smiles in ground_truth]

        if not self._canonicalized:
            predictions = [
                [canonicalize_smiles(smiles) for smiles in smiles_list]
                for smiles_list in predictions
            ]
            ground_truth = [
                [canonicalize_smiles(smiles) for smiles in smiles_list]
                for smiles_list in ground_truth
            ]

        is_in_set = []
        for i in range(len(ground_truth)):
            ground_truth_i = ground_truth[i]
            predictions_i = predictions[i]

            nr_tgt_smi = 0
            count = 0
            for tgt_smi in ground_truth_i:
                nr_tgt_smi += 1
                if tgt_smi in predictions_i:
                    count += 1

            coverage = count / nr_tgt_smi

            is_in_set.append(coverage)

        return is_in_set

    def _score_predictions(
        self, predictions: PredictionType, ground_truth: TargetType
    ) -> Dict[str, float]:
        is_in_set = self._is_in_set(predictions, ground_truth)

        score = np.mean(is_in_set)

        return {self.scorer_name: score}
