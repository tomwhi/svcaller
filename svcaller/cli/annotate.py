import json

import click

import pandas as pd

from svcaller.effect.consequence import predict_effects
from svcaller.calling.events import read_sv_gtf, SvType


@click.command()
@click.argument('svs', type=click.Path(exists=True))
@click.option('--ts-regions', type=click.Path(exists=True), required=True)
@click.option('--ar-regions', type=click.Path(exists=True), required=True)
@click.option('--fusion-regions', type=click.Path(exists=True), required=True)
@click.option('--effects-filename', type=str, default="sv_effects.json", required=False)
@click.pass_context
def predict_effects_cmd(ctx, svs, ts_regions, ar_regions, fusion_regions, effects_filename):
    with open(effects_filename, 'w') as effects_file, open(svs) as svs_file, \
        open(ts_regions) as ts_regions_file, open(ar_regions) as ar_regions_file, \
        open(fusion_regions) as fusion_regions_file:
        effects = predict_effects(svs_file, ts_regions_file, ar_regions_file, fusion_regions_file)
        json.dump(effects, effects_file)


@click.command()
@click.argument('svs-bed', type=str)
@click.option('--del-gtf', type=click.Path(exists=True), required=True)
@click.option('--dup-gtf', type=click.Path(exists=True), required=True)
@click.option('--inv-gtf', type=click.Path(exists=True), required=True)
@click.option('--tra-gtf', type=click.Path(exists=True), required=True)
@click.pass_context
def gtf_to_bed_cmd(ctx, svs_bed, del_gtf, dup_gtf, inv_gtf, tra_gtf):
    with open(svs_bed, 'w') as svs_bed_file, open(del_gtf) as del_gtf_file, \
        open(dup_gtf) as dup_gtf_file, open(inv_gtf) as inv_gtf_file, \
        open(tra_gtf) as tra_gtf_file:
        gtf_files_and_event_types = \
            zip([del_gtf_file, dup_gtf_file, inv_gtf_file, tra_gtf_file],
                [SvType.DEL, SvType.DUP, SvType.INV, SvType.TRA])
        sv_bed_tables = [read_sv_gtf(gtf_file, event_type)
                         for (gtf_file, event_type) in gtf_files_and_event_types]
        svs_bed_table = pd.concat(sv_bed_tables)
        json.dump(svs_bed_file, svs_bed_table)