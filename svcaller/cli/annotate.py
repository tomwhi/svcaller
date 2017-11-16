import json

import click

from svcaller.effect.consequence import predict_effects


@click.command()
@click.argument('svs', type=click.Path(exists=True))
@click.option('--ts-regions', type=click.Path(exists=True), required=True)
@click.option('--ar-regions', type=click.Path(exists=True), required=True)
@click.option('--fusion-regions', type=click.Path(exists=True), required=True)
@click.option('--effects-filename', type=str, default="sv_effects.json", required=False)
@click.pass_context
def predict_effects_cmd(ctx, svs, ts_regions, ar_regions, fusion_regions, effects_filename):
    with open(effects_filename, 'w') as effects_file:
        effects = predict_effects(svs, ts_regions, ar_regions, fusion_regions)
        json.dump(effects, effects_file)
