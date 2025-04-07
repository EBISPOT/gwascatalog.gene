from typing import Annotated

from pydantic import BaseModel, PositiveFloat, Field, PositiveInt, AliasChoices


GeneName = Annotated[
    str,
    Field(description="Gene name", validation_alias=AliasChoices("Name", "gene_name")),
]

Chromosome = Annotated[str, Field(description="Chromosome name")]

Position = Annotated[
    PositiveInt,
    Field(
        description="Base pair position",
        validation_alias=AliasChoices("base_pair_location"),
    ),
]

pValue = Annotated[
    PositiveFloat,
    Field(
        le=1,
        description="p-value of GWAS association",
        validation_alias=AliasChoices("p_value"),
    ),
]


class GeneModel(BaseModel):
    """A row in a gene-based sumstat file"""
    name: GeneName
    chromosome: Chromosome
    position: Position
    pValue: pValue
