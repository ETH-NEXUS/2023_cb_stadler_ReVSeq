import {TableDataSampleCounts, TableDataMetadata} from 'src/models/tables'
import {SampleCount, Metadata} from 'src/models/core'

export const createSampleCountsData = (sampleCounts: SampleCount[]): TableDataSampleCounts[] => {
  console.log('sampleCounts', sampleCounts)
  return sampleCounts.map(sampleCount => {
    return {
      pseudonymized_id: sampleCount.sample
        ? sampleCount.sample.pseudonymized_id
        : sampleCount.pseudonymized_id,
      substrain: sampleCount.substrain.name,
      strain: sampleCount.substrain.strain.name,
      aligned: sampleCount.aligned,
      length: sampleCount.length,
      rpkm: sampleCount.rpkm,
      rpkm_proportions: sampleCount.rpkm_proportions,
      outlier: sampleCount.outlier,
      coverage_threshold: sampleCount.coverage_threshold,
      coverage: sampleCount.coverage,
      coverage_status: sampleCount.coverage_status,
      readnum_threshold: sampleCount.readnum_threshold,
      readnum_status: sampleCount.readnum_status,
      percentile_threshold: sampleCount.percentile_threshold,
      tax_id: sampleCount.tax_id,
      scientific_name: sampleCount.scientific_name,
      DP: sampleCount.DP,
      consensus_fraction_n: sampleCount.consensus_fraction_n,
      consensus_number_n: sampleCount.consensus_number_n,
      consensus: sampleCount.consensus,
      consensus_cds: sampleCount.consensus_cds,
      mean_coverage_non_N_positions: sampleCount.mean_coverage_non_N_positions,



    }
  })
}

export const createMetadataRows = (metadata: Metadata[]): TableDataMetadata[] => {
  return metadata.map(m => {
    return {
      pseudonymized_id: m.sample.pseudonymized_id,
      prescriber: m.prescriber,
      well: m.well.location,
      order_date: m.order_date,
      ent_date: m.ent_date,
      treatment_type: m.treatment_type,
    }
  })
}
