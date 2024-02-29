export interface SampleCountsTableColumns {
  name: string
  label: string
  field: string | ((row: any) => any)
  required?: boolean
  align?: 'left' | 'right' | 'center'
  sortable?: boolean
  sort?: (a: any, b: any, rowA: any, rowB: any) => number
  sortOrder?: 'ad' | 'da'
  format?: (val: any, row: any) => any
  style?: string | ((row: any) => string)
  classes?: string | ((row: any) => string) | undefined
}
export const columnsMetadata: SampleCountsTableColumns[] = [
  {
    name: 'pseudonymized_id',
    label: 'Sample ID',
    field: 'pseudonymized_id',
    align: 'center',
    sortable: true,
  },
  {name: 'well', label: 'Well', field: 'well', align: 'left', sortable: true},
  {name: 'prescriber', label: 'Prescriber Kanton', field: 'prescriber', align: 'left', sortable: true},
  {name: 'order_date', label: 'Order date', field: 'order_date', align: 'left', sortable: true},
  {name: 'ent_date', label: 'Ent date', field: 'ent_date', align: 'left', sortable: true},
  {name: 'treatment_type', label: 'Treatment type', field: 'treatment_type', align: 'left', sortable: true},
]

export const columnsSampleCount: SampleCountsTableColumns[] = [
  {
    name: 'pseudonymized_id',
    label: 'Sample ID',
    field: 'pseudonymized_id',
    align: 'center',
    sortable: true,
  },
  {
    name: 'substrain',
    required: true,
    label: 'Substrain',
    align: 'left',
    sortable: true,
    field: 'substrain',
  },
  {
    name: 'strain',
    required: true,
    label: 'Strain',
    align: 'center',
    sortable: true,
    field: 'strain',
  },
  {
    name: 'aligned',
    required: true,
    label: 'Aligned',
    align: 'center',
    sortable: true,
    field: 'aligned',
  },
  {
    name: 'length',
    required: true,
    label: 'Length',
    align: 'center',
    sortable: true,
    field: 'length',
  },
  {
    name: 'rpkm',
    required: true,
    label: 'RPKM',
    align: 'center',
    sortable: true,
    field: 'rpkm',
  },
  {
    name: 'rpkm_proportions',
    required: true,
    label: 'RPKM proportions',
    align: 'center',
    sortable: true,
    field: 'rpkm_proportions',
  },
  {
    name: 'normcounts',
    required: true,
    label: 'Normalized counts',
    align: 'center',
    sortable: true,
    field: 'normcounts',
  },
  {
    name: 'outlier',
    required: true,
    label: 'Outlier',
    align: 'center',
    sortable: true,
    field: 'outlier',
  },
  {
    name: 'coverage',
    required: true,
    label: 'Coverage',
    align: 'center',
    sortable: true,
    field: 'coverage',
  },
  {
    name: 'coverage_threshold',
    required: true,
    label: 'Coverage threshold',
    align: 'center',
    sortable: true,
    field: 'coverage_threshold',

  },
  {
    name: 'coverage_status',
    required: true,
    label: 'Coverage status',
    align: 'center',
    sortable: true,
    field: 'coverage_status',

  },
  {
    name: 'readnum_threshold',
    required: true,
    label: 'Readnum threshold',
    align: 'center',
    sortable: true,
    field: 'readnum_threshold',

  },
  {
    name: 'readnum_status',
    required: true,
    label: 'Readnum status',
    align: 'center',
    sortable: true,
    field: 'readnum_status',

  },
  {
    name: 'percentile_threshold',
    required: true,
    label: 'Percentile threshold',
    align: 'center',
    sortable: true,
    field: 'percentile_threshold',

  },

  {
    name: 'DP20',
    required: true,
    label: 'DP20',
    align: 'center',
    sortable: true,
    field: 'DP20',
  }
  ,

  {
    name: 'DP1',
    required: true,
    label: 'DP1',
    align: 'center',
    sortable: true,
    field: 'DP1',
  }
  ,

  {
    name: 'DP2',
    required: true,
    label: 'DP2',
    align: 'center',
    sortable: true,
    field: 'DP2',
  }

]
